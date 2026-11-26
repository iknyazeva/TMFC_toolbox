function tmfc_LSS_estimate(SPM,fmriData,defaults,isRWLS,outNii)

%--------------------------------------------------------------------------
% Estimate LSS GLM and save the current-trial beta image
%--------------------------------------------------------------------------
% This is an adaptation of the original SPM25 GLM estimation procedure.
% Functional data are reused from a run-level cache or streamed in chunks.
% Temporal non-sphericity is estimated independently for each LSS trial,
% using SPM ReML or rWLS ReML. Only the first beta image is saved.
%
% For the default non-rWLS F-contrast, ESS is computed as rssB - ResSS,
% where rssB is calculated once per session and cached.
% For rWLS and custom F-contrasts, ESS is computed directly.
% ReML estimation, whitening, and beta estimation are performed
% independently for each trial.
%
% FORMAT tmfc_LSS_estimate(SPM,fmriData,defaults,isRWLS,outNii)
% SPM      - LSS GLM structure for the current trial
% fmriData - cached functional data or information for streaming from disk
% defaults - TMFC defaults (including the ReML threshold ufp)
% isRWLS   - true for robust weighted least squares estimation
% outNii   - output path for the current-trial beta image
%__________________________________________________________________________
%
% References:
%
% Statistical Parametric Maps in Functional Imaging: A General Linear
% Approach. Friston KJ, Holmes AP, Worsley KJ, Poline JB, Frith CD,
% Frackowiak RSJ. (1995) Human Brain Mapping 2:189-210.
%
% Analysis of fMRI Time-Series Revisited - Again. Worsley KJ, Friston KJ.
% (1995) NeuroImage 2:173-181.
%
% Diedrichsen & Shadmehr (2006) Detecting and Adjusting for artifacts in
%       fMRI time series data
%__________________________________________________________________________
% Original SPM code:
% Karl Friston & Guillaume Flandin
% Copyright (C) 1994-2022 Wellcome Centre for Human Neuroimaging
%
% rWLS: Diedrichsen & Shadmehr (2006)
%
%
% TMFC modifications:
% Copyright (C) 2026 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

%-Get design
%--------------------------------------------------------------------------
xX = SPM.xX;
xVi = SPM.xVi;
nScan = size(xX.X,1);

%-Check consistency with cached functional data
%--------------------------------------------------------------------------
if nScan ~= fmriData.nscan
    error('tmfc_LSS:ScanMismatch','Design and data have different scan counts.');
end
if ~isequaln(xX.K,fmriData.K)
    error('tmfc_LSS:FilterMismatch','Trial filter differs from the cached run-level filter.');
end
if numel(xX.name) < 2 || isempty(strfind(xX.name{1},'Current_trial')) || ...
        isempty(strfind(xX.name{2},'Other_trials'))
    error('tmfc_LSS:Design','The first two regressors must be Current_trial / Other_trials.');
end

%-Estimate temporal non-sphericity (ReML) for the current trial
%==========================================================================
selectedMask = cell(numel(fmriData.blocks),1);
if ~isfield(xVi,'V')

    %-First pass with identity whitening
    %----------------------------------------------------------------------
    xKXs = spm_sp('Set',spm_filter(xX.K,xX.X));
    xKXs.X = full(xKXs.X);

    %-Effects-of-interest F-contrast
    %----------------------------------------------------------------------
    if isfield(xVi,'Fcontrast')
        fc = spm_FcUtil('Set','User-specified contrast','F','c',xVi.Fcontrast,xKXs);
    else
        fc = spm_FcUtil('Set','effects of interest','F','iX0',[xX.iB xX.iG],xKXs);
    end

    % Use the cached ESS baseline for the default non-rWLS F-contrast
    % when the baseline design matches the current trial
    iX0 = [xX.iB xX.iG];
    useRssB = ~isRWLS && ~isfield(xVi,'Fcontrast') && ...
        ~isempty(iX0) && ~isempty(fc(1).c) && ...
        isfield(fmriData,'essBaseline') && ...
        isequaln(fmriData.essBaseline.iX0,iX0) && ...
        isequaln(fmriData.essBaseline.X0,xX.X(:,iX0));

    A = [];
    if isempty(fc(1).c)
        trMV = 1;
    else
        X1o = spm_FcUtil('X1o',fc(1),xKXs);
        trMV = spm_SpUtil('trMV',X1o);
        if ~useRssB
            % Prepare the projection matrix for direct ESS computation
            pKX0 = spm_sp('x-',xKXs);
            Hsqr = spm_FcUtil('Hsqr',fc(1),xKXs);
            A = Hsqr*pKX0;
        end
    end
    trRV = spm_SpUtil('trRV',xKXs);
    UF = spm_invFcdf(1-defaults.ufp,[trMV,trRV]);
    if isRWLS, UF = 0; end  % Disable F-thresholding for rWLS

    %-Select voxels and accumulate spatial covariance
    %----------------------------------------------------------------------
    Cy = zeros(nScan,nScan);
    nSelected = 0;
    for b=1:numel(fmriData.blocks)
        Y = tmfc_LSS_chunk_data(fmriData,b);
        if isempty(Y), continue; end
        if fmriData.cached
            KY = fmriData.blocks(b).KY;
        else
            KY = spm_filter(xX.K,Y);
        end
        ResSS = sum(spm_sp('r',xKXs,KY).^2,1);
        if useRssB
            ESS = fmriData.essBaseline.rssB{b} - ResSS;
            j = ESS/trMV > UF*ResSS/trRV;
        elseif isempty(A)
            j = true(1,size(Y,2));
        else
            ESS = sum((A*KY).^2,1);
            j = ESS/trMV > UF*ResSS/trRV;
        end
        if ~isRWLS
            % Restrict the final beta mask to voxels selected for ReML
            selectedMask{b} = j;
        end
        if any(j)
            weights = sqrt(trRV./ResSS(j));
            % Avoid unnecessary indexing when all voxels are selected
            if all(j)
                Y = bsxfun(@times,Y,weights);
            else
                Y = bsxfun(@times,Y(:,j),weights);
            end
            Cy = Cy + Y*Y';
            nSelected = nSelected + nnz(j);
        end
    end
    if nSelected == 0
        error('tmfc_LSS:ReMLMask','No voxels selected for trial-wise ReML.');
    end
    Cy = Cy/nSelected;

    %-Estimate covariance separately for each session block
    %----------------------------------------------------------------------
    if isstruct(xX.K)
        V = sparse(nScan,nScan);
        for s=1:numel(xX.K)
            rows = xX.K(s).row;
            Q = {};
            for j=1:numel(xVi.Vi)
                Qj = xVi.Vi{j}(rows,rows);
                if nnz(Qj)
                    Q{end+1} = Qj;
                end
            end
            Xp = xX.X(rows,:);
            if isfield(xX.K(s),'X0')
                Xp = [Xp xX.K(s).X0];
            end
            if isRWLS
                [Vp,~] = tmfc_spm_rwls_reml(Cy(rows,rows),Xp,Q);
            else
                [Vp,~] = spm_reml(Cy(rows,rows),Xp,Q);
            end
            V(rows,rows) = V(rows,rows)+Vp;
        end
    else
        if isRWLS
            [V,~] = tmfc_spm_rwls_reml(Cy,xX.X,xVi.Vi);
        else
            [V,~] = spm_reml(Cy,xX.X,xVi.Vi);
        end
    end
    V = V*nScan/trace(V);
else
    V = xVi.V;
end

%-Get weight/whitening matrix: W*W' = inv(V)
%--------------------------------------------------------------------------
W = spm_sqrtm(spm_inv(V));
W = W.*(abs(W)>1e-6);

%-Construct only the first beta row of the GLS estimator
%--------------------------------------------------------------------------
xKXs = spm_sp('Set',spm_filter(xX.K,W*xX.X));
xKXs.X = full(xKXs.X);
pKX = spm_sp('x-',xKXs);
row = (W'*spm_filter(xX.K,pKX(1,:)'))';

%-Initialise the current-trial beta image (beta_0001)
%--------------------------------------------------------------------------
Vbeta = struct('fname',outNii, ...
    'dim',fmriData.dim,'dt',[spm_type('float32') spm_platform('bigend')], ...
    'mat',fmriData.mat,'pinfo',[1 0 0]', ...
    'descrip',sprintf('spm_spm:beta (%04d) - %s',1,xX.name{1}));
Vbeta = spm_data_hdr_write(Vbeta);

%-Estimate beta coefficients and write the image in chunks
%==========================================================================
for b=1:numel(fmriData.blocks)
    ind = fmriData.blocks(b).indices;
    ix = fmriData.blocks(b).range;
    c = NaN(numel(ix),1);
    if ~isempty(ind)
        Y = tmfc_LSS_chunk_data(fmriData,b);
        beta = row*Y;
        if ~isRWLS && ~isempty(selectedMask{b})
            valid = selectedMask{b};
            ind = ind(valid);
            beta = beta(valid);
        end
        c(ind-ix(1)+1) = beta(:);
    end
    Vbeta = spm_data_write(Vbeta,c,ix);
end
end

%--------------------------------------------------------------------------
% Read a chunk of functional data from cache or disk
%--------------------------------------------------------------------------
function Y = tmfc_LSS_chunk_data(fmriData,b)
    if fmriData.cached
        Y = fmriData.blocks(b).Y;
        return;
    end
    ix = fmriData.blocks(b).indices;
    if isempty(ix), Y = zeros(fmriData.nscan,0); return; end
    Y = zeros(fmriData.nscan,numel(ix));
    for j=1:fmriData.nscan
        Y(j,:) = spm_data_read(fmriData.V(j),ix);
    end
end
