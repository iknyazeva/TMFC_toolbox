function Cy = tmfc_LSS_pooled_cy(SPM,D,R,UFp)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Compute the run-pooled sample covariance for LSS using the nuisance-only
% model (session regressors and constants). Replaces the per-trial sample
% covariance accumulation in SPM's spm_est_non_sphericity. Unlike main,
% this covariance is shared; ReML and whitening remain trial-specific.
%
% FORMAT Cy = tmfc_LSS_pooled_cy(SPM,D,R,UFp)
% SPM has xX.X set to R.N; D and R hold the run data and filtered data.
% UFp is the SPM voxel-pooling threshold.
%
% =========================================================================
% Derived from spm_est_non_sphericity.m (SPM25) and tmfc_LSS.m
% License: GPL-3.0-or-later

xX = SPM.xX;

%-Design-level quantities under i.i.d. (spm_est_non_sphericity.m:59-92)
xKXs   = spm_sp('Set',spm_filter(xX.K,xX.X));
xKXs.X = full(xKXs.X);
iX0    = [];
xCon   = spm_FcUtil('Set','effects of interest','F','iX0',iX0,xKXs);
if isempty(xCon(1).c)
    error('tmfc_LSS_pooled_cy: run design without effects of interest is not supported.');
end
X1o  = spm_FcUtil('X1o',xCon(1),xKXs);
trMV = spm_SpUtil('trMV',X1o);
trRV = spm_SpUtil('trRV',xKXs);
UF   = spm_invFcdf(1 - UFp,[trMV,trRV]);

%-Pooled variance estimation (:128-198)
Cy    = 0;
nPool = 0;
for i = 1:numel(D.Y)
    Y = D.Y{i};
    if isempty(Y), continue, end

    KY    = R.KY{i};
    ResSS = sum(spm_sp('r',xKXs,KY).^2,1);          % spm_est_non_sphericity.m:175-180
    rssB  = sum(KY.^2,1);
    ESS   = rssB - ResSS;

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
end
