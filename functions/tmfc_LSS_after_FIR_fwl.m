function Fw = tmfc_LSS_after_FIR_fwl(KXt,sN)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% LSS-after-FIR helper: Frisch-Waugh-Lovell basis of the trial
% columns.
%
% The LSS-after-FIR model is y = [x_t N] [b_t; b_N] + e, where x_t holds
% the columns of Current_trial and Other_trials (q = 2) and N holds
% everything else (here only the session constant(s): ones(nScan,1), or
% blkdiag(ones(n_1,1),...,ones(n_S,1)) for concatenated runs, because the
% after-FIR models carry no confounds). With the filtered design (K leaves
% the data unchanged here because hpf = Inf, but the call is kept to mirror
% main), FWL gives the same least-squares fit as the full model:
%
%   rN   = (I - P_KN)*KY                        residual of the data on N
%   Z    = (I - P_KN)*K*x_t                     residual of the trial columns
%   b_t  = pinv(Z)*rN
%
% P_KN and rN do not depend on the trial and are computed once per run
% (tmfc_LSS_after_FIR_fwl_batch); only Z (n x q) is computed per trial.
%
% Rank decision: Z is scaled column-wise by the norm of the original filtered
% column K*x_t, so the singular values of Zs = Z./norm are the fractions of each
% column that are not explained by N, and the tolerance max(size)*eps(max(sv))
% (the rule of spm_sp) is invariant to the scale of the columns. A column that
% lies in span(KN) or in the span of the other trial column gives a singular
% value of the order of eps.
%
% FORMAT Fw = tmfc_LSS_after_FIR_fwl(KXt,sN)
%
%   KXt - [n x q] filtered trial columns
%   sN  - spm_sp space of the filtered nuisance columns
%
%   Fw.Q   - orthonormal basis of span(Z) [n x rk]
%   Fw.Vs  - V(:,1:rk)./s(1:rk)': Zs pseudo-inverse is Fw.Vs*Fw.Q'
%   Fw.sc  - column norms of KXt
%   Fw.rk  - rank of Z
%   Fw.sv  - singular values of Zs
%   Fw.deficient - rk < q
%
% =========================================================================
% License: GPL-3.0-or-later

q  = size(KXt,2);
Z  = spm_sp('r',sN,KXt);
sc = sqrt(sum(KXt.^2,1));
sc(sc == 0) = 1;
[U,S,V] = svd(Z./sc,'econ');
s   = diag(S);
tol = max(size(Z))*eps(max(s));
rk  = sum(s > tol);

Fw.Q   = U(:,1:rk);
Fw.Vs  = V(:,1:rk)./s(1:rk)';
Fw.sc  = sc;
Fw.rk  = rk;
Fw.sv  = s;
Fw.deficient = rk < q;
end
