function [beta1,W] = tmfc_LSS_pooled_estimate(SPM,D,Cy,R,h0,keep,W)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Pooled LSS estimator (lss_method = 'pooled', "LSS with run-pooled
% covariance"): estimate the first beta (Current_trial) of one LSS GLM with
% an approximate separation of the per-trial ReML. Design, data handling,
% whitening construction (W = spm_sqrtm(spm_inv(V)), |W| <= 1e-6 -> 0),
% pseudo-inverse (spm_sp SVD tolerance) and the beta row solve are those of
% the exact estimator (tmfc_LSS_exact_estimate); the difference to main is
% the model change: spm_reml(Cy, X_k, Vi) uses the run-level Cy of
% tmfc_LSS_pooled_cy instead of a per-trial Cy_k, so no per-trial first
% pass is run. Per-trial cost: ReML + pKX + n*v row solve (the n^2*v Cy
% accumulation is gone). The estimate is approximate by design; the
% accuracy against main is measured by tests/lss/run_lss_tests.m.
%
% i.i.d. parents (xVi.V given): the parent W is used, as in the exact
% estimator; nothing is approximated.
%
% FORMAT beta1 = tmfc_LSS_pooled_estimate(SPM,D,Cy,R)
%
%   SPM - design-only SPM structure of one trial (tmfc_LSS_design)
%   D   - run data (tmfc_LSS_load_data)
%   Cy  - run-pooled covariance from tmfc_LSS_pooled_cy
%   R   - run-level part from tmfc_LSS_exact_prepare (R.W0 for i.i.d.)
%
%   beta1 - {nChunk x 1} first beta in the in-mask voxels of each chunk
%
% =========================================================================
% Derived from tmfc_LSS.m (Copyright (C) 2026 Ruslan Masharipov) and SPM25
% License: GPL-3.0-or-later

if nargin < 5, h0 = []; end
if nargin < 6, keep = []; end
if nargin < 7, W = []; end
xX  = SPM.xX;
xVi = SPM.xVi;

if ~isempty(W)
    % Reuse the last estimated whitening matrix within this run.
elseif isfield(xVi,'V')
    % i.i.d.: V is given (parent W), nothing is approximated
    W = R.W0;
else
    % ReML on the run-pooled Cy (spm_est_non_sphericity.m:245-285)
    [V,~] = reml_block(xX,xVi,Cy,h0,keep);
    V     = V*size(xX.X,1)/trace(V);
    W     = spm_sqrtm(spm_inv(V));
    W     = W.*(abs(W) > 1e-6);
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
function [V,h] = reml_block(xX,xVi,Cy,h0,keep)
% spm_reml per session block (spm_est_non_sphericity.m:245-285)

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
