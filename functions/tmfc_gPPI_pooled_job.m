function job = tmfc_gPPI_pooled_job(SPM,sess_num)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Whole-brain load job for the gPPI pooled-Cy path (gppi_method='pooled',
% analysis 2): an fmri_spec-shaped job whose sess.scans is the scan list of
% ALL analyzed runs in design order, i.e. the data set of the seed designs
% (the same list tmfc_gPPI/tmfc_gPPI_FIR put into
% base_spec.sess(kSess).scans / the per-seed spec batches). Feeding the
% whole concatenated list to tmfc_LSS_load_data in ONE call resolves the
% cross-run question of OPTIONS_DESIGN R1: spm_est_non_sphericity.m:56
% (nScan = size(xX.X,1), the concatenated timeline) and :133-158 pool over
% the concatenated data with one joint per-chunk voxel mask, so the shared
% Cy must be estimated over all runs jointly (per-run loads would produce
% per-run masks and chunk grids that do not align across runs).
%
% The caller must invoke this WHILE the PARENT SPM is loaded (later
% per-seed loads overwrite the SPM variable) and must reject concatenated
% sessions beforehand (SPM_concat == 1 lists every image in every session,
% which would duplicate scans).
%
% FORMAT job = tmfc_gPPI_pooled_job(SPM,sess_num)
%
%   job     - fmri_spec-shaped job for tmfc_LSS_load_data (sess.scans,
%             global, mthresh, mask; the other spec fields are unused by
%             the loader)
%   SPM     - parent SPM of the subject (xY.VY, Sess, xGX, xM)
%   sess_num- analyzed session (run) numbers, in design order
%
% License: GPL-3.0-or-later. Derived from the TMFC toolbox (R. Masharipov)
% and SPM25 spm_est_non_sphericity.m.

job = struct();
job.global  = SPM.xGX.iGXcalc;
job.mthresh = SPM.xM.gMT;
try
    job.mask = {SPM.xM.VM.fname};
catch
    job.mask = {''};
end

scans = cell(0,1);
for kSess = 1:numel(sess_num)
    row = SPM.Sess(sess_num(kSess)).row;
    for image = 1:numel(row)
        VYi = SPM.xY.VY(row(image));
        scans{end + 1,1} = [VYi.fname ',' num2str(VYi.n(1))]; %#ok<AGROW>
    end
end
job.sess.scans = scans;
end
