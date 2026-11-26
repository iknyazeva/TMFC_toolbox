function pKX = tmfc_gPPI_estimate(SPM,beta_idx,beta_paths,tmfc)

%--------------------------------------------------------------------------
% Estimate gPPI GLM and save selected PPI beta images
%--------------------------------------------------------------------------
% This is an adaptation of the original SPM25 function spm_spm.
% It preserves SPM's GLM estimation procedure but saves only selected
% beta images, without residual images, smoothness estimates, or SPM.mat.
%
% FORMAT pKX = tmfc_gPPI_estimate(SPM,beta_idx,beta_paths,tmfc)
% SPM        - gPPI GLM structure
% beta_idx   - column indices of selected PPI regressors
% beta_paths - output paths for selected beta images
% tmfc       - TMFC configuration structure
% pKX        - pseudoinverse of the filtered and whitened design matrix
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
%__________________________________________________________________________
% Original SPM code:
% Karl Friston & Guillaume Flandin
% Copyright (C) 1994-2022 Wellcome Centre for Human Neuroimaging
%
%
% TMFC modifications:
% Copyright (C) 2026 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

if numel(beta_idx) ~= numel(beta_paths)
    error('tmfc_gPPI:BetaSelection','beta_idx and beta_paths must have same length.');
end
if numel(unique(beta_idx)) ~= numel(beta_idx)
    error('tmfc_gPPI:BetaSelection','Duplicate beta indices are not supported.');
end
VY = SPM.xY.VY;
for i = 1:numel(VY)
    if ~spm_existfile(VY(i).fname)
        error('File not found: %s',VY(i).fname);
    end
    if ~spm_mesh_detect(VY)
        % Backward compatibility: propagate scaling (see spm_fmri_spm_ui.m)
        VY(i).private.dat.scl_slope = VY(i).pinfo(1);
        VY(i).private.dat.scl_inter = VY(i).pinfo(2);
    end
end
spm_check_orientations(VY);
M = VY(1).mat;
DIM = VY(1).dim;
YNaNrep = spm_type(VY(1).dt(1),'nanrep');
if spm_mesh_detect(VY)
    file_ext = '.gii';
    g = VY(1).private;
    metadata = g.private.metadata;
    name = {metadata.name};
    if any(ismember(name,'SurfaceID'))
        metadata = metadata(ismember(name,'SurfaceID'));
        metadata = {metadata.name, metadata.value};
    elseif isfield(g,'faces') && ~isempty(g.faces)
        metadata = {'SurfaceID',VY(1).fname};
    else
        error('SurfaceID not found in GIfTI''s metadata.');
    end
    if isempty(spm_file(metadata{2},'path'))
        metadata{2} = fullfile(spm_file(VY(1).fname,'path'),metadata{2});
    end
else
    file_ext = spm_file_ext;
    metadata = {};
end

%-Get design
%--------------------------------------------------------------------------
xX             = SPM.xX;
[nScan, nBeta] = size(xX.X);

%-Get masking settings
%--------------------------------------------------------------------------
if isfield(SPM,'xM')
    xM         = SPM.xM;
else
    xM         = -Inf(nScan,1);
end
if ~isstruct(xM)
    xM         = struct(...
                    'T',  [],...
                    'TH', xM,...
                    'I',  0,...
                    'VM', {[]},...
                    'xs', struct('Masking','analysis threshold'));
end

mask           = true(DIM);

%-Check confounds (xX.K)
%--------------------------------------------------------------------------
if ~isfield(xX,'K')
    xX.K       = 1;
end

%-Get non-sphericity (xVi), otherwise assume i.i.d.
%--------------------------------------------------------------------------
if isfield(SPM,'xVi')
    xVi        = SPM.xVi;
else
    xVi        = struct('form', 'i.i.d.',...
                        'V',    speye(nScan,nScan));
end

%-Evoke ReML for hyperparameter estimation
%--------------------------------------------------------------------------
if ~isfield(xVi,'V')
    SPM.xY.VY  = VY;
    SPM.xM     = xM;
    SPM.xX.K   = xX.K;
    [xVi, am]  = tmfc_spm_est_non_sphericity(SPM,tmfc);
    mask       = mask & am;
    spm('FnBanner',mfilename);
end

%-Get weight/whitening matrix:  W*W' = inv(V)
%--------------------------------------------------------------------------
if isfield(xX,'W')
    W          = xX.W;
else
    W          = spm_sqrtm(spm_inv(xVi.V));
    W          = W.*(abs(W) > 1e-6);
    xX.W       = sparse(W);
end

%-Design space and projector matrix [pseudoinverse] for WLS
%--------------------------------------------------------------------------
xX.xKXs        = spm_sp('Set',spm_filter(xX.K,W*xX.X));    % KWX
xX.xKXs.X      = full(xX.xKXs.X);
xX.pKX         = spm_sp('x-',xX.xKXs);                     % Projector
erdf           = spm_SpUtil('trRV',xX.xKXs);               % error df

%-Use non-sphericity xVi.V to compute [effective] degrees of freedom
%--------------------------------------------------------------------------
xX.V           = spm_filter(xX.K,spm_filter(xX.K,W*xVi.V*W')'); % KWVW'K'
[trRV, trRVRV] = spm_SpUtil('trRV',xX.xKXs,xX.V);          % trRV (for X)
xX.trRV        = trRV;                                     % <R'*y'*y*R>
xX.trRVRV      = trRVRV;                                   %-Satterthwaite
xX.erdf        = trRV^2/trRVRV;                            % approximation
xX.Bcov        = xX.pKX*xX.V*xX.pKX';                      % Cov(beta)

%-Initialise only the requested beta files
%--------------------------------------------------------------------------
if any(beta_idx < 1) || any(beta_idx > nBeta) || any(beta_idx ~= fix(beta_idx))
    error('tmfc_gPPI:BetaSelection','Beta column index is outside model design.');
end
Vbeta = repmat(struct( ...
    'fname', [], 'dim', DIM, ...
    'dt', [spm_type('float32') spm_platform('bigend')], ...
    'mat', M, 'pinfo', [1 0 0]', 'descrip', 'spm_spm:beta', ...
    metadata{:}), 1, numel(beta_idx));
for i = 1:numel(beta_idx)
    Vbeta(i).fname = beta_paths{i};
    Vbeta(i).descrip = sprintf('spm_spm:beta (%04d) - %s',beta_idx(i),xX.name{beta_idx(i)});
end
if ~isempty(Vbeta)
    Vbeta = spm_data_hdr_write(Vbeta);
end

%-Get explicit mask(s)
%==========================================================================
for i = 1:numel(xM.VM)
    if ~(isfield(SPM,'xVol') && isfield(SPM.xVol,'G'))
        %-Assume it fits entirely in memory
        C = spm_bsplinc(xM.VM(i), [0 0 0 0 0 0]');
        v = true(DIM);
        [x1,x2] = ndgrid(1:DIM(1),1:DIM(2));
        for x3 = 1:DIM(3)
            M2  = inv(M\xM.VM(i).mat);
            y1 = M2(1,1)*x1+M2(1,2)*x2+(M2(1,3)*x3+M2(1,4));
            y2 = M2(2,1)*x1+M2(2,2)*x2+(M2(2,3)*x3+M2(2,4));
            y3 = M2(3,1)*x1+M2(3,2)*x2+(M2(3,3)*x3+M2(3,4));
            v(:,:,x3) = spm_bsplins(C, y1,y2,y3, [0 0 0 0 0 0]') > 0;
        end
        mask = mask & v;
        clear C v x1 x2 x3 M2 y1 y2 y3
    else
        if spm_mesh_detect(xM.VM(i))
            v = full(xM.VM(i).private.cdata) > 0;
        else
            v = spm_mesh_project(gifti(SPM.xVol.G), xM.VM(i)) > 0;
        end
        mask = mask & v(:);
        clear v
    end
end

%-Split data into chunks
%==========================================================================
chunksize = floor(tmfc.defaults.maxmem / 8 / nScan);
nbchunks  = ceil(prod(DIM) / chunksize);
chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),prod(DIM)+1);

for i=1:nbchunks
    chunk = chunks(i):chunks(i+1)-1;
    
    %-Get data & construct analysis mask
    %======================================================================
    Y     = zeros(nScan,numel(chunk));
    cmask = mask(chunk);
    for j=1:nScan
        if ~any(cmask), break, end                 %-Break if empty mask
        
        Y(j,cmask) = spm_data_read(VY(j),chunk(cmask));%-Read chunk of data
        
        cmask(cmask) = Y(j,cmask) > xM.TH(j);      %-Threshold (& NaN) mask
        if xM.I && ~YNaNrep && xM.TH(j) < 0        %-Use implicit mask
            cmask(cmask) = abs(Y(j,cmask)) > eps;
        end
    end
    cmask(cmask) = any(diff(Y(:,cmask),1));        %-Mask constant data
    Y            = Y(:,cmask);                     %-Data within mask
    
    %-Whiten/Weight data and remove filter confounds
    %======================================================================
    KWY          = spm_filter(xX.K,W*Y);
    
    %-Weighted Least Squares estimation
    %======================================================================
    beta         = xX.pKX*KWY;                     %-Parameter estimates
    % Residuals and residual-based smoothness statistics are not needed
    % for first-level PPI beta images used by TMFC group analysis.

    %-Write only PPI beta images
    %----------------------------------------------------------------------
    mask(chunk) = cmask;
    c = NaN(numel(chunk),1);
    for j = 1:numel(beta_idx)
        c(cmask) = beta(beta_idx(j),:);
        Vbeta(j) = spm_data_write(Vbeta(j),c,chunk);
    end
end

if ~any(mask(:))
    error('Please check your data: There are no inmask voxels.');
end
pKX = xX.pKX;
end
