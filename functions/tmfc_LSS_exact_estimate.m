function beta1 = tmfc_LSS_exact_estimate(SPM,D,R,UFp)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Exact LSS estimator (lss_method = 'exact'): estimate the first beta
% (Current_trial) of one LSS GLM with the SAME per-trial ReML as tmfc_LSS of
% main (branch accelerated_lss); only the linear algebra on the data differs:
%
%   - first pass (W = I): ResSS_k = sum(spm_sp('r',xKXs,KY).^2) with K*Y
%     hoisted out of the trial loop (tmfc_LSS_exact_prepare), bit-identical
%     to the per-trial computation of main; ESS of the 'effects of interest'
%     F-contrast = rssB - ResSS_k (rssB is run-level), so no beta and no
%     Hsqr*beta on the data;
%   - pool voxels, weights sqrt(trRV/ResSS) and Cy = Y*diag(w)*Y' as in
%     spm_est_non_sphericity.m;
%   - ReML, W = spm_sqrtm(spm_inv(V)), |W| <= 1e-6 -> 0: as in main;
%   - second pass: beta_1 = g*Y with g = pKX(1,:)*K*W (one row, n x 1
%     products instead of W*Y, spm_filter and pKX*KWY on the data). pKX is
%     the spm_sp pseudo-inverse of K*W*X as in spm_spm.m:438-440, so rank
%     deficiency and the SVD tolerance are those of tmfc_LSS of main.
%
% The difference to main is the order of equivalent floating-point
% operations on the data: the outputs are numerically identical to main
% within float32 output precision (to be checked by tests/lss/run_lss_tests.m).
%
% FORMAT beta1 = tmfc_LSS_exact_estimate(SPM,D,R,UFp)
%
%   SPM  - design-only SPM structure of one trial (tmfc_LSS_design)
%   D    - run data (tmfc_LSS_load_data)
%   R    - run-level part from tmfc_LSS_exact_prepare
%   UFp  - voxel threshold p-value as in spm_est_non_sphericity (1 in tmfc_LSS)
%
%   beta1 - {nChunk x 1} first beta in the in-mask voxels of each chunk
%
% =========================================================================
% Derived from tmfc_LSS.m (Copyright (C) 2026 Ruslan Masharipov) and SPM25
% License: GPL-3.0-or-later

xX  = SPM.xX;
xVi = SPM.xVi;

if ~isequal(xX.X(:,3:end),R.N)
    error('tmfc_LSS_exact_estimate: nuisance part of the design differs from the run-level part.');
end

%-First pass (spm_spm.m:417-424)
if ~isfield(xVi,'V')
    xVi = reml_pass(xX,xVi,D,R,UFp);
end

%-Whitening matrix (spm_spm.m:431-433)
if ~isempty(R.W0)
    W = R.W0;
else
    W = spm_sqrtm(spm_inv(xVi.V));
    W = W.*(abs(W) > 1e-6);
end

%-Design space and pseudo-inverse of K*W*X (spm_spm.m:438-440)
xKXs   = spm_sp('Set',spm_filter(xX.K,W*xX.X));
xKXs.X = full(xKXs.X);
pKX    = spm_sp('x-',xKXs);

%-beta_1 = pKX(1,:)*K*(W*Y) = (pKX(1,:)*K*W)*Y
g     = (W'*spm_filter(xX.K,pKX(1,:)'))';
beta1 = cell(numel(D.Y),1);
for i = 1:numel(D.Y)
    beta1{i} = g*D.Y{i};
end
end

%==========================================================================
function xVi = reml_pass(xX,xVi,D,R,UFp)
% First pass of spm_est_non_sphericity.m (W = I) with the run-level K*Y of
% tmfc_LSS_exact_prepare: identical operations and order to the per-trial
% computation of main (ResSS from spm_sp, ESS = R.rssB - ResSS)

nScan = size(xX.X,1);

%-Design-level quantities under i.i.d. (:59-92), n x m objects
xKXs   = spm_sp('Set',spm_filter(xX.K,xX.X));
xKXs.X = full(xKXs.X);
xCon   = spm_FcUtil('Set','effects of interest','F','iX0',R.iX0,xKXs);
if isempty(xCon(1).c)
    error('tmfc_LSS_exact_estimate: design without effects of interest is not supported.');
end
X1o  = spm_FcUtil('X1o',xCon(1),xKXs);
trMV = spm_SpUtil('trMV',X1o);
trRV = spm_SpUtil('trRV',xKXs);
xVi.UFp = UFp;
UF   = spm_invFcdf(1 - UFp,[trMV,trRV]);

%-Pooled variance estimation (:128-198)
Cy    = 0;
nPool = 0;
for i = 1:numel(D.Y)
    Y = D.Y{i};
    if isempty(Y), continue, end

    ResSS = sum(spm_sp('r',xKXs,R.KY{i}).^2,1);          % spm_est_non_sphericity.m:175-180
    ESS   = R.rssB{i} - ResSS;

    j     = ESS/trMV > UF*ResSS/trRV;
    q     = nnz(j);
    nPool = nPool + q;
    if q
        w = sqrt(trRV./ResSS(j));
        if q == size(Y,2)
            Y = Y.*w;
        else
            Y = Y(:,j).*w;
        end
        Cy = Cy + Y*Y';
    end
end

if ~nPool
    error('Please check your data: There are no significant voxels.');
end
Cy = Cy / nPool;

%-ReML estimation (:245-285)
if isstruct(xX.K)
    m = length(xVi.Vi);
    h = zeros(m,1);
    V = sparse(nScan,nScan);
    for i = 1:length(xX.K)
        q  = xX.K(i).row;
        p  = [];
        Qp = {};
        for j = 1:m
            if nnz(xVi.Vi{j}(q,q))
                Qp{end + 1} = xVi.Vi{j}(q,q); %#ok<AGROW>
                p           = [p j]; %#ok<AGROW>
            end
        end
        Xp = xX.X(q,:);
        try %#ok<TRYNC>
            Xp = [Xp xX.K(i).X0]; %#ok<AGROW>
        end
        [Vp,hp] = spm_reml(Cy(q,q),Xp,Qp);
        V(q,q)  = V(q,q) + Vp; %#ok<SPRIX>
        h(p)    = hp;
    end
else
    [V,h] = spm_reml(Cy,xX.X,xVi.Vi);
end

V     = V*nScan/trace(V);
xVi.h = h;
xVi.V = V;
end
