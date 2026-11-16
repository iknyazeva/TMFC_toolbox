function D = tmfc_LSS_after_FIR_data(job)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% LSS-after-FIR helper: read the residual data (Res_*.nii of the FIR
% task regression) of one run once and keep the in-mask time series in
% memory. The after-FIR model uses the following data handling:
%
%   - NO global scaling of the FIR residuals: gSF = 1 for both
%     global 'None' and 'Scaling' (FIR residual scaling correction). spm_global is not called, the
%     VY pinfo / scl_slope / scl_inter of the Res_*.nii headers are used
%     as they are on disk. (The original routine re-applied the
%     grand-mean scaling of spm_fmri_spm_ui.m:299-343 to the residuals;
%     the fixed reference removes that after fmri_spec, so the data read
%     here is bit-identical to what spm_spm reads in the fixed main.)
%   - threshold mask TH = -Inf(nScan,1) for every scan: as in main, which
%     forces SPM.xM.TH = -Inf(size(SPM.xM.TH)) after fmri_spec (main uses
%     mthresh = -Inf; g.*gSF*(-Inf) is NOT used because it would give NaN
%     or +Inf for g <= 0). Y > -Inf keeps every voxel except NaN voxels.
%   - explicit mask (SPM.swd mask.nii, spm_bsplinc/spm_bsplins > 0),
%     constant-voxel removal (any(diff(Y,1))) and maxmem chunking are the
%     SPM operations (spm_est_non_sphericity.m:101-166, spm_spm.m:527-585)
%     as in the LSS input loader.
%
% Everything here is independent of the trial, whereas main recomputes it
% for every trial inside spm_jobman('run', ...).
%
% FORMAT D = tmfc_LSS_after_FIR_data(job)
%
%   job - matlabbatch{1}.spm.stats.fmri_spec of any trial of the run
%         (the fields used here are the same for all trials of the run;
%         job.mthresh must be -Inf)
%
%   D.Y     - {nChunk x 1} in-mask data [nScan x nVoxelsInChunk], unscaled
%   D.cmask - {nChunk x 1} in-mask index within the chunk
%   D.chunk - {nChunk x 1} linear voxel indices of the chunk
%   D.DIM, D.M - image dimensions and voxel-to-world matrix of the first scan
%   D.nIn   - number of in-mask voxels
%   D.gSF   - ones(nScan,1) (no scaling; kept for diagnostics/parity)
%   D.TH    - -Inf(nScan,1) (see above)
%
% Chunks are the same as in spm_spm (stats.maxmem), so voxel sums and
% writes run in the same order as in main.
%
% =========================================================================
% Derived from tmfc_LSS.m (Copyright (C) 2026 Ruslan Masharipov), from
% SPM25.
% License: GPL-3.0-or-later

% Scans as in spm_run_fmri_spec.m:129
scans = job.sess.scans;
P     = strvcat(scans{:});
nScan = numel(scans);
if nScan == 1
    error('Not enough scans in session 1.');
end

% after-FIR models always use mthresh = -Inf
if ~(isinf(job.mthresh) && job.mthresh < 0)
    error('tmfc_LSS_after_FIR_data: job.mthresh must be -Inf for LSS-after-FIR.');
end

VY = spm_data_hdr_read(P);
spm_check_orientations(VY);

% No global scaling of the FIR residuals:
% gSF = 1 for both 'None' and 'Scaling'; the header scale factors of the
% Res_*.nii files are used unchanged, exactly as restored in the fixed main.
gSF = ones(nScan,1);

% Masking threshold: -Inf for every scan (main sets SPM.xM.TH = -Inf after
% fmri_spec; Y > -Inf removes only NaN voxels of the residuals)
TH = -Inf(nScan,1);

M   = VY(1).mat;
DIM = VY(1).dim;

% Explicit mask (spm_est_non_sphericity.m:101-126)
mask = true(DIM);
if ~isempty(job.mask{1})
    VM = spm_data_hdr_read(job.mask{1});
    for i = 1:numel(VM)
        C = spm_bsplinc(VM(i), [0 0 0 0 0 0]');
        v = true(DIM);
        [x1,x2] = ndgrid(1:DIM(1),1:DIM(2));
        for x3 = 1:DIM(3)
            M2  = inv(M\VM(i).mat);
            y1 = M2(1,1)*x1+M2(1,2)*x2+(M2(1,3)*x3+M2(1,4));
            y2 = M2(2,1)*x1+M2(2,2)*x2+(M2(2,3)*x3+M2(2,4));
            y3 = M2(3,1)*x1+M2(3,2)*x2+(M2(3,3)*x3+M2(3,4));
            v(:,:,x3) = spm_bsplins(C, y1,y2,y3, [0 0 0 0 0 0]') > 0;
        end
        mask = mask & v;
    end
end

% Chunks (spm_est_non_sphericity.m:133-135)
chunksize = floor(spm_get_defaults('stats.maxmem') / 8 / nScan);
nbchunks  = ceil(prod(DIM) / chunksize);
chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),prod(DIM)+1);

D.Y     = cell(nbchunks,1);
D.cmask = cell(nbchunks,1);
D.chunk = cell(nbchunks,1);
D.nIn   = 0;
for i = 1:nbchunks
    chunk = chunks(i):chunks(i+1)-1;

    % Get data and construct analysis mask (spm_est_non_sphericity.m:150-166)
    Y     = zeros(nScan,numel(chunk));
    cmask = mask(chunk);
    for j = 1:nScan
        if ~any(cmask), break, end
        Y(j,cmask) = spm_data_read(VY(j),chunk(cmask));
        cmask(cmask) = Y(j,cmask) > TH(j);
    end
    cmask(cmask) = any(diff(Y(:,cmask),1));
    Y = Y(:,cmask);

    D.Y{i}     = Y;
    D.cmask{i} = cmask;
    D.chunk{i} = chunk;
    D.nIn      = D.nIn + nnz(cmask);
end

D.DIM = DIM;
D.M   = M;
D.gSF = gSF;
D.TH  = TH;

if D.nIn == 0
    error('Please check your data: There are no inmask voxels.');
end
end
