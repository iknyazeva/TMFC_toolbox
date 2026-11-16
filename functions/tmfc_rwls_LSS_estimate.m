function beta1 = tmfc_rwls_LSS_estimate(SPM,D,UFp,KY,Cy)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% rWLS LSS estimator: estimate the first beta
% (Current_trial) of one LSS GLM in memory with the SAME rWLS algorithm as
% the rWLS branch of tmfc_LSS of main (tmfc_rwls_LSS + tmfc_spm_rwls_spm,
% Diedrichsen & Shadmehr 2006). Used by the accelerated tmfc_LSS
% (the accelerated implementation) for rWLS parents (SPM.xVi.form = 'wls'); the
% output is compared to main with a relative-error gate
% (tests/lss/run_lss_tests.m):
%
%   1. first pass at W = I on the in-memory run data: OLS beta, residuals,
%      ResSS, voxel pool and weights sqrt(trRV/ResSS), pooled spatially
%      whitened covariance Cy (tmfc_spm_rwls_est_non_sphericity.m:57-232;
%      the wls branch forces xVi.UFp = inf, UF = 0, so selection uses
%      positive F-contrast energy, including the original edge cases);
%   2. one tmfc_spm_rwls_reml per session block of xX.K on the per-scan
%      covariance components (tmfc_spm_rwls_est_non_sphericity.m:245-285,
%      the wls branch); V is normalised to trace nScan;
%   3. W = spm_sqrtm(spm_inv(V)) with |W| <= 1e-6 set to zero
%      (tmfc_spm_rwls_spm.m / spm_spm.m:431-433);
%   4. K*W*X pseudo-inverse by spm_sp; apply its first row to K*W*Y
%      (tmfc_spm_rwls_spm.m:438-440, 589-593), reassociating the product.
%
% Not computed (not needed for beta_0001): residuals/ResStats of the second
% pass, ResMS, RPV, smoothness, xX.V, trRV/trRVRV/erdf, further betas.
%
% SPM must carry the rWLS covariance structure exactly as the rWLS branch of
% tmfc_LSS of main sets it: xVi.form = 'wls' and one component Vi{j} per scan
% (sparse selector with a single 1 at (j,j)) - see the caller.
%
% FORMAT beta1 = tmfc_rwls_LSS_estimate(SPM,D,UFp,KY,Cy)
%
%   SPM  - design-only SPM structure (tmfc_LSS_design) with xVi set to
%          the rWLS per-scan components of the run
%   D    - run data (tmfc_LSS_load_data)
%   UFp  - pool threshold as in tmfc_LSS (1); in the wls branch UF is
%          forced to 0, so the value does not change the result
%   KY   - optional run-level filtered data
%   Cy   - optional approximate run-pooled covariance; ReML is per trial
%
%   beta1 - {nChunk x 1} first beta in the in-mask voxels of each chunk
%
% =========================================================================
% Derived from tmfc_LSS.m / denoising/tmfc_spm_rwls_spm.m /
% denoising/tmfc_spm_rwls_est_non_sphericity.m (rWLS toolbox v4.1)
% License: GPL-3.0-or-later

xX   = SPM.xX;
xVi  = SPM.xVi;

if ~strcmp(xVi.form,'wls')
    error('tmfc_rwls_LSS_estimate: SPM.xVi.form must be ''wls'' (per-scan components).');
end

%-First pass and ReML (tmfc_spm_rwls_est_non_sphericity.m, W = I)
if nargin < 4, KY = []; end
if nargin < 5, Cy = []; end
xVi = rwls_pass(xX,xVi,D,UFp,KY,Cy);

%-Get weight/whitening matrix: W*W' = inv(V) (tmfc_spm_rwls_spm.m:431-433)
W = spm_sqrtm(spm_inv(xVi.V));
W = W.*(abs(W) > 1e-6);

%-Design space and projector matrix for WLS (tmfc_spm_rwls_spm.m:438-440)
xX.xKXs   = spm_sp('Set',spm_filter(xX.K,W*xX.X));
xX.xKXs.X = full(xX.xKXs.X);
xX.pKX    = spm_sp('x-',xX.xKXs);

%-Only the first beta row is needed. Apply the same operator to the row
% before multiplying by the data, as in tmfc_LSS_exact_estimate.
g = (W'*spm_filter(xX.K,xX.pKX(1,:)'))';
beta1 = cell(numel(D.Y),1);
for i = 1:numel(D.Y)
    beta1{i} = g*D.Y{i};
end
end

%==========================================================================
function xVi = rwls_pass(xX,xVi,D,UFp,KY,Cy)
% First pass of tmfc_spm_rwls_est_non_sphericity on the in-memory run data,
% or an optional run-pooled Cy, followed by per-block rWLS ReML.

nScan = size(xX.X,1);
if isempty(Cy)
xX.W  = speye(nScan,nScan);

%-Compute Hsqr and F-threshold under i.i.d. (:57-92)
xX.xKXs   = spm_sp('Set',spm_filter(xX.K,xX.W*xX.X));
xX.xKXs.X = full(xX.xKXs.X);
xX.pKX    = spm_sp('x-',xX.xKXs);

iX0  = [xX.iB xX.iG];
xCon = spm_FcUtil('Set','effects of interest','F','iX0',iX0,xX.xKXs);
if ~isempty(xCon(1).c)
    X1o  = spm_FcUtil('X1o', xCon(1),xX.xKXs);
    Hsqr = spm_FcUtil('Hsqr',xCon(1),xX.xKXs);
    trMV = spm_SpUtil('trMV',X1o);
else
    trMV = 1;
    Hsqr = Inf;
end
trRV = spm_SpUtil('trRV',xX.xKXs);

%-Threshold for voxels entering the estimates (:85-102); the wls branch
% forces xVi.UFp = inf, UF = 0; retain the original selection expression.
xVi.UFp = UFp;
UF      = spm_invFcdf(1 - UFp,[trMV,trRV]); %#ok<NASGU>
xVi.UFp = inf;
UF      = 0;

%-Pooled variance estimation (:128-232)
Cy    = 0;
nPool = 0;
for i = 1:numel(D.Y)
    Y = D.Y{i};
    if isempty(Y), continue, end

    if isempty(KY)
        KWY = spm_filter(xX.K,xX.W*Y);
    else
        KWY = KY{i};
    end
    beta  = xX.pKX*KWY;
    res   = spm_sp('r',xX.xKXs,KWY);
    ResSS = sum(res.^2);
    clear res

    j     = sum((Hsqr*beta).^2,1)/trMV > UF*ResSS/trRV;
    q     = nnz(j);
    nPool = nPool + q;
    if q
        q  = spdiags(sqrt(trRV./ResSS(j)'),0,q,q);
        Y  = Y(:,j)*q;
        Cy = Cy + Y*Y';
    end
end

if ~nPool
    error('Please check your data: There are no significant voxels.');
end
Cy = Cy / nPool;
end

%-ReML estimation per session block (:245-285, wls branch ->
% tmfc_spm_rwls_reml instead of spm_reml)
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
        [Vp,hp] = tmfc_spm_rwls_reml(Cy(q,q),Xp,Qp);
        V(q,q)  = V(q,q) + Vp; %#ok<SPRIX>
        h(p)    = hp;
    end
else
    [V,h] = tmfc_spm_rwls_reml(Cy,xX.X,xVi.Vi);
end

%-Normalize non-sphericity and save hyperparameters (:287-290)
V     = V*nScan/trace(V);
xVi.h = h;
xVi.V = V;
xVi.Cy = Cy;
end
