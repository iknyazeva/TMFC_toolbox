function beta1 = tmfc_rwls_LSS(SPM,D,UFp,KY,Cy)
% Estimate one rWLS LSS trial in memory, replacing the local tmfc_rwls_LSS
% and tmfc_spm_rwls_spm job in main. Uses rWLS ReML per trial; the first
% beta row is evaluated with the same operator in a different arithmetic order.
% Copyright (C) 2026 Ruslan Masharipov. GPL-3.0-or-later.
if isfield(SPM.xVi,'V')
    SPM.xVi = rmfield(SPM.xVi,'V');
end
SPM.xVi.form = 'wls';
nScan = size(SPM.xX.X,1);
SPM.xVi.Vi = cell(nScan,1);
for jScan = 1:nScan
    SPM.xVi.Vi{jScan} = sparse(jScan,jScan,1,nScan,nScan);
end
if nargin < 4, KY = []; end
if nargin < 5, Cy = []; end
beta1 = tmfc_rwls_LSS_estimate(SPM,D,UFp,KY,Cy);
end
