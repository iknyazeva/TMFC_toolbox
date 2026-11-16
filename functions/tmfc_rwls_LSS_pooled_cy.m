function Cy = tmfc_rwls_LSS_pooled_cy(SPM,D,KY)
% Run-pooled rWLS covariance from the nuisance-only design. The voxel
% selection and variance scaling follow tmfc_spm_rwls_est_non_sphericity;
% this Cy is approximate for trial designs, whose ReML remains separate.
% Derived from tmfc_spm_rwls_est_non_sphericity.m. GPL-3.0-or-later.
xX = SPM.xX;
xX.X = xX.X(:,3:end);
xX.iB = xX.iB - 2; xX.iG = xX.iG - 2;
s = spm_sp('Set',spm_filter(xX.K,xX.X)); s.X = full(s.X);
p = spm_sp('x-',s);
c = spm_FcUtil('Set','effects of interest','F','iX0',[xX.iB xX.iG],s);
if ~isempty(c(1).c)
    X1o = spm_FcUtil('X1o',c(1),s);
    Hsqr = spm_FcUtil('Hsqr',c(1),s);
    trMV = spm_SpUtil('trMV',X1o);
else
    trMV = 1; Hsqr = Inf;
end
trRV = spm_SpUtil('trRV',s);
Cy = 0; nPool = 0;
for i = 1:numel(D.Y)
    Y = D.Y{i}; if isempty(Y), continue, end
    beta = p*KY{i};
    res = spm_sp('r',s,KY{i});
    ResSS = sum(res.^2);
    j = sum((Hsqr*beta).^2,1)/trMV > 0*ResSS/trRV;
    q = nnz(j); nPool = nPool + q;
    if q
        w = spdiags(sqrt(trRV./ResSS(j)'),0,q,q);
        Y = Y(:,j)*w;
        Cy = Cy + Y*Y';
    end
end
if ~nPool, error('Please check your data: There are no significant voxels.'); end
Cy = Cy/nPool;
end
