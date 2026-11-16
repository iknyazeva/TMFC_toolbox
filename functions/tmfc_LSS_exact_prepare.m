function R = tmfc_LSS_exact_prepare(SPM,D)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Exact LSS helper (lss_method = 'exact'): the trial-independent part of the
% first (W = I) ReML pass, computed once per run by the accelerated tmfc_LSS
% (branch accelerated_lss). Replaces the per-trial recomputation of these
% quantities in tmfc_LSS of main; identical operations and order, so the
% exact estimates are numerically identical to main within float32 output
% precision.
%
% FORMAT R = tmfc_LSS_exact_prepare(SPM,D)
%
%   SPM - design-only SPM of any trial of the run (tmfc_LSS_design)
%   D   - run data (tmfc_LSS_load_data)
%
%   R.N    - nuisance columns of the design (all except Current/Other)
%   R.iX0  - columns of xX.iB and xX.iG (constants, DCT)
%   R.KY   - {nChunk} K*Y (hoisted out of the trial loop; the filter does not
%            depend on the trial)
%   R.rssB - {nChunk} residual sum of squares of K*Y on K*X0, X0 = columns of
%            xX.iB and xX.iG; ESS of the F-contrast 'effects of interest' is
%            rssB - ResSS
%   R.W0   - whitening matrix for i.i.d. parents (xVi.V given), else []
%
% =========================================================================
% Derived from tmfc_LSS.m (Copyright (C) 2026 Ruslan Masharipov) and SPM25
% License: GPL-3.0-or-later

xX = SPM.xX;
if ~contains(xX.name{1},'Current_trial') || ~contains(xX.name{2},'Other_trials')
    error('tmfc_LSS_exact_prepare: columns 1 and 2 of the design must be Current_trial and Other_trials.');
end

R.N    = xX.X(:,3:end);
R.iX0  = [xX.iB xX.iG];
sB     = spm_sp('Set',spm_filter(xX.K,xX.X(:,R.iX0)));

nC = numel(D.Y);
R.KY = cell(nC,1); R.rssB = cell(nC,1);
for i = 1:nC
    KY       = spm_filter(xX.K,D.Y{i});
    R.KY{i}  = KY;
    R.rssB{i} = sum(spm_sp('r',sB,KY).^2,1);
end

% i.i.d.: V is given, W is the same for all trials (spm_spm.m:431-433)
R.W0 = [];
if isfield(SPM.xVi,'V')
    W    = spm_sqrtm(spm_inv(SPM.xVi.V));
    W    = W.*(abs(W) > 1e-6);
    R.W0 = sparse(W);
end
end
