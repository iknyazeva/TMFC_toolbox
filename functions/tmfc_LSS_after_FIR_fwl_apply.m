function [beta1,w] = tmfc_LSS_after_FIR_fwl_apply(Fw,rN)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% LSS-after-FIR helper: trial weights and the first (Current_trial)
% beta from the FWL basis (tmfc_LSS_after_FIR_fwl) and the residual of the
% data on the nuisance part.
% LSS-after-FIR needs only the first beta row (no ResSS, no ReML), so the
% weights of the first column are returned explicitly.
%
% FORMAT [beta1,w] = tmfc_LSS_after_FIR_fwl_apply(Fw,rN)
%
%   Fw   - basis from tmfc_LSS_after_FIR_fwl
%   rN   - [n x v] residual of the data on the nuisance part, spm_sp('r',sN,KY)
%
%   beta1 - [1 x v] Current_trial coefficient: w*rN
%   w     - [1 x n] row of the FWL pseudo-inverse for the first column,
%           w = (Fw.Vs(1,:)/Fw.sc(1))*Fw.Q'
%
% beta1 = (Fw.Vs(1,:)*Fw.Q'*rN)/Fw.sc(1) is the first row of pinv(Z)*rN in
% the column-normalised basis, the FWL equivalent of pKX(1,:)*K*Y in spm_spm.
%
% =========================================================================
% License: GPL-3.0-or-later

w = (Fw.Vs(1,:)./Fw.sc(1))*Fw.Q';
if nargin > 1 && ~isempty(rN)
    beta1 = w*rN;
else
    beta1 = [];
end
end
