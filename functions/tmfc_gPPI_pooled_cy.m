function [Cy,nPool] = tmfc_gPPI_pooled_cy(SPM,D,UFp)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Run/seed-pooled sample covariance for gPPI / gPPI-FIR analysis 2
% (gppi_method='pooled'). Adapted from tmfc_LSS_pooled_cy.m (single-run
% trial design, branch accelerated_lss) and SPM25
% spm_est_non_sphericity.m:66-77,128-198: the nuisance-only model (the
% seed design minus the PPI_* and Seed_* regressor columns; FIR condition
% regressors stay) replaces the full seed design, and the voxel pool and
% Cy are estimated ONCE per subject over the concatenated timeline of all
% analyzed runs and shared across seeds. ReML and whitening stay per seed
% (tmfc_gPPI_pooled_estimate). The estimate is approximate by design
% (OPTIONS_DESIGN section 3 / Revision R1).
%
% Like spm_est_non_sphericity, Cy accumulates Y*Y' of the SELECTED,
% variance-scaled RAW (unfiltered) time series, including cross-run
% blocks; the downstream reml_block consumes only the per-session diagonal
% blocks Cy(q,q) (spm_est_non_sphericity.m:292), while the pooled voxel
% count normalizes the whole matrix jointly (spm_est_non_sphericity.m:197).
%
% FORMAT [Cy,nPool] = tmfc_gPPI_pooled_cy(SPM,D,UFp)
%
%   SPM   - STRIPPED seed SPM of the first ROI (xX.W / xVi.V removed by
%           the caller); xX spans the concatenated timeline of all runs
%   D     - whole-brain data of ALL runs from ONE tmfc_LSS_load_data call
%           over the concatenated scan list (tmfc_gPPI_pooled_job):
%           D.Y{i} is [nScan_total x nVoxelsInChunk] (spatial chunks)
%   UFp   - voxel-pooling threshold p (SPM stats.fmri.ufp)
%   Cy    - [nScan_total x nScan_total] pooled sample covariance
%   nPool - number of pooled (selected) voxels
%
% License: GPL-3.0-or-later. Derived from the TMFC toolbox (R. Masharipov)
% and SPM25 spm_est_non_sphericity.m.

xX = SPM.xX;

%-Nuisance-only design: drop the seed-specific PPI/PHYS columns
% (names are ['PPI_' PPI.PPI.name] and ['Seed_' ROIname]; see the
% regress construction in tmfc_gPPI.m / tmfc_gPPI_FIR.m). SPM25 fmri_spec
% prefixes EVERY xX.name entry with the session tag 'Sn(k) ' (multi-run:
% 'Sn(1)'..'Sn(k)'; xX.name is a cellstr), so the tag is stripped before
% the prefix match -- a raw strncmp(xX.name,...) never matched (s16d).
nm = regexprep(xX.name,'^Sn\(\d+\)\s+','');
isPPI = strncmp(nm,'PPI_',4) | strncmp(nm,'Seed_',5);
if ~any(isPPI)
    error('tmfc_gPPI_pooled_cy:noPpiColumns', ...
        'no PPI_/Seed_ columns found in the seed design');
end
ikeep = find(~isPPI);
xX.X    = xX.X(:,ikeep);
xX.name = xX.name(ikeep);
[allB,xX.iB] = ismember(xX.iB,ikeep); %#ok<ASGLU>  % constants/DCT are never PPI/Seed
[allG,xX.iG] = ismember(xX.iG,ikeep); %#ok<ASGLU>
if ~all(allB) || ~all(allG)
    error('tmfc_gPPI_pooled_cy:columnMismatch', ...
        'iB/iG columns must survive the PPI/Seed column drop');
end

%-Design-level quantities under i.i.d. (spm_est_non_sphericity.m:66-77)
xX.xKXs   = spm_sp('Set',spm_filter(xX.K,xX.X));
xX.xKXs.X = full(xX.xKXs.X);
xX.pKX    = spm_sp('x-',xX.xKXs);
xCon      = spm_FcUtil('Set','effects of interest','F','iX0',[xX.iB xX.iG],xX.xKXs);
if isempty(xCon(1).c)
    error('tmfc_gPPI_pooled_cy:noEffects', ...
        'run design without effects of interest is not supported');
end
X1o  = spm_FcUtil('X1o',xCon(1),xX.xKXs);
Hsqr = spm_FcUtil('Hsqr',xCon(1),xX.xKXs);
trMV = spm_SpUtil('trMV',X1o);
trRV = spm_SpUtil('trRV',xX.xKXs);
UF   = spm_invFcdf(1 - UFp,[trMV,trRV]);

%-Pooled variance estimation over spatial chunks
% (spm_est_non_sphericity.m:133-166 data pass and :168-198 selection)
Cy    = 0;
nPool = 0;
for i = 1:numel(D.Y)
    Y = D.Y{i};
    if isempty(Y), continue, end

    KY    = spm_filter(xX.K,Y);
    beta  = xX.pKX*KY;
    ResSS = sum(spm_sp('r',xX.xKXs,KY).^2,1);

    j = sum((Hsqr*beta).^2,1)/trMV > UF*ResSS/trRV;
    q = nnz(j);
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
Cy = Cy/nPool;
end
