function D = tmfc_LSS_load_data(job)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Read the data of one run once and keep the in-mask time series in memory
% (run-level inputs for the accelerated tmfc_LSS).
%
% Everything here is independent of the trial, whereas tmfc_LSS of main
% recomputes it for every trial (spm_fmri_spm_ui.m + the data pass of
% spm_est_non_sphericity.m and spm_spm.m); this function does it once per
% run:
%
%   - file handles, global signal, grand-mean scaling, threshold mask
%     (spm_fmri_spm_ui.m:285-364)
%   - explicit mask (spm_est_non_sphericity.m:101-126)
%   - threshold / NaN / constant-voxel mask and chunking
%     (spm_est_non_sphericity.m:131-166, spm_spm.m:553-585)
%
% FORMAT D = tmfc_LSS_load_data(job)
%
%   job - matlabbatch{1}.spm.stats.fmri_spec of any trial of the run
%         (the fields used here are the same for all trials of the run)
%
%   D.Y     - {nChunk x 1} in-mask data [nScan x nVoxelsInChunk], scaled
%   D.cmask - {nChunk x 1} in-mask index within the chunk
%   D.chunk - {nChunk x 1} linear voxel indices of the chunk
%   D.DIM, D.M - image dimensions and voxel-to-world matrix of the first scan
%   D.nIn   - number of in-mask voxels
%   D.g, D.gSF, D.TH - global signal, scaling factor and mask threshold per scan
%
% Chunks are the same as in spm_spm (stats.maxmem), so sums over voxels run
% in the same order as in tmfc_LSS of main; the outputs are numerically
% identical to main within float32 output precision.
%
% =========================================================================
% Derived from tmfc_LSS.m (Copyright (C) 2026 Ruslan Masharipov) and SPM25
% License: GPL-3.0-or-later

% Scans as in spm_run_fmri_spec.m:129
scans = job.sess.scans;
P     = strvcat(scans{:}); %#ok<DSTRVCT>
nScan = numel(scans);
if nScan == 1
    error('Not enough scans in session 1.');
end

VY = spm_data_hdr_read(P);
spm_check_orientations(VY);

% Global signal and grand-mean scaling (spm_fmri_spm_ui.m:299-343)
GM  = 100;
g   = zeros(nScan,1);
for i = 1:nScan
    g(i) = spm_global(VY(i));
end
gSF = GM./g;
if strcmpi(job.global,'none')
    gSF(1:nScan) = GM./mean(g);
end
for i = 1:nScan
    VY(i).pinfo(1:2,:) = VY(i).pinfo(1:2,:) * gSF(i);
    VY(i).private.dat.scl_slope = VY(i).private.dat.scl_slope * gSF(i);
    VY(i).private.dat.scl_inter = VY(i).private.dat.scl_inter * gSF(i);
end

% Masking threshold (spm_fmri_spm_ui.m:364)
TH = g.*gSF*job.mthresh;

% Propagate scaling (spm_spm.m:309-318)
for i = 1:nScan
    VY(i).private.dat.scl_slope = VY(i).pinfo(1);
    VY(i).private.dat.scl_inter = VY(i).pinfo(2);
end

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
D.g   = g;
D.gSF = gSF;
D.TH  = TH;

if D.nIn == 0
    error('Please check your data: There are no inmask voxels.');
end
end
