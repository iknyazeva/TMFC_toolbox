function [beta,W,h] = tmfc_gPPI_pooled_estimate(SPM,Y,Cy,h0,keep)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Pooled gPPI estimator (gppi_method='pooled', analysis 2): estimate the
% ROI-to-ROI betas of ONE seed with an approximate separation of the
% per-seed ReML (OPTIONS_DESIGN Revision R2). Structure follows
% tmfc_LSS_pooled_estimate.m:33-61 (branch accelerated_lss): ReML on the
% run-pooled Cy per K-block, V*nScan/trace(V) normalization, whitening
% W = spm_sqrtm(spm_inv(V)) with |W| <= 1e-6 set to zero, pseudo-inverse
% pKX of K*W*X. Unlike the LSS helper (which solves only the first beta
% row), the FULL beta = pKX*Y is returned because the a2 caller consumes
% every column of the preloaded ROI stack Y. The whitening matrix W between
% seeds is NOT reused (an extra approximation outside this slice); with
% warm start only the shared Cy and the h0 cache carry over (R4).
%
% FORMAT [beta,W,h] = tmfc_gPPI_pooled_estimate(SPM,Y,Cy,h0,keep)
%
%   SPM  - STRIPPED seed SPM (xX.W / xVi.V removed by the caller, exact
%          path parity); must NOT carry xX.W or xVi.V
%   Y    - preloaded ROI data stack of the caller ([nScan_total x nROI])
%   Cy   - run-pooled covariance from tmfc_gPPI_pooled_cy
%   h0   - optional warm-start hyperparameters (empty = default SPM
%          initialisation); expected from the first estimated seed (R4)
%   keep - optional logical mask of retained Vi components (h0 ~= 0)
%
%   beta - [nBeta x nROI] beta matrix (pKX*Y, all rows)
%   W    - whitening matrix of this seed
%   h    - hyperparameters of this seed (h0 cache for the next seeds)
%
% License: GPL-3.0-or-later. Derived from the TMFC toolbox (R. Masharipov)
% and SPM25 spm_est_non_sphericity.m.

if nargin < 4, h0 = []; end
if nargin < 5, keep = []; end

xX  = SPM.xX;
xVi = SPM.xVi;
if isfield(xX,'W') || isfield(xVi,'V')
    error('tmfc_gPPI_pooled_estimate:notStripped', ...
        'caller must strip xX.W/xVi.V first (exact path parity)');
end

%-ReML on the run-pooled Cy (spm_est_non_sphericity.m:245-285)
[V,h] = reml_block(xX,xVi,Cy,h0,keep);
V     = V*size(xX.X,1)/trace(V);
W     = spm_sqrtm(spm_inv(V));
W     = W.*(abs(W) > 1e-6);

%-Design space and pseudo-inverse of K*W*X (spm_spm.m:438-440)
xKXs   = spm_sp('Set',spm_filter(xX.K,W*xX.X));
xKXs.X = full(xKXs.X);
pKX    = spm_sp('x-',xKXs);

%-beta = pKX*Y (FULL stack; the a2 caller whitens every ROI column)
beta = pKX*Y;
end

%==========================================================================
function [V,h] = reml_block(xX,xVi,Cy,h0,keep)
% Local copy of the reml_block subfunction of tmfc_LSS_pooled_estimate.m
% (:66-105, branch accelerated_lss; the LSS helper itself is NOT
% modified): spm_reml per session block, optional warm start via
% tmfc_spm_reml_warm (a warm-start copy of SPM25 spm_reml: "only the
% initial h may differ").

nScan = size(xX.X,1);
m     = length(xVi.Vi);
if isstruct(xX.K)
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
        if ~isempty(keep)
            selected = keep(p); p = p(selected); Qp = Qp(selected);
        end
        Xp = xX.X(q,:);
        try %#ok<TRYNC>
            Xp = [Xp xX.K(i).X0]; %#ok<AGROW>
        end
        if isempty(h0)
            [Vp,hp] = spm_reml(Cy(q,q),Xp,Qp);
        else
            [Vp,hp] = tmfc_spm_reml_warm(Cy(q,q),Xp,Qp,[],[],[],[],h0(p));
        end
        V(q,q)  = V(q,q) + Vp; %#ok<SPRIX>
        h(p)    = hp;
    end
else
    p = 1:m; Qp = xVi.Vi;
    if ~isempty(keep), p = p(keep); Qp = Qp(keep); end
    if isempty(h0)
        [V,hp] = spm_reml(Cy,xX.X,Qp);
    else
        [V,hp] = tmfc_spm_reml_warm(Cy,xX.X,Qp,[],[],[],[],h0(p));
    end
    h = zeros(m,1); h(p) = hp;
end
end
