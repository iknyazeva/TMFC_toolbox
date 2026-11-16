function results = run_lss_after_fir_synth_tests(out_dir,sets,spm_dir,inputs_root)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Compare the in-memory and batched LSS-after-FIR estimators with the fixed
% reference on 10 synthetic input sets:
%   <inputs_root>/
%       {02,02v,12,15,17}/Sub_0{1,2}/parent/SPM.mat
% For every set the script runs the FIR task regression once (tmfc_FIR,
% window 24 s, 12 bins, serial) into
%   <out_dir>/<design>/<Sub_0k>/fir/FIR_regression/<Sub_0k>/Res_*.nii
% and then runs the reference and both estimators on the same inputs,
% serially with default BLAS threads and no parfor.
%
% Comparison: ONLY the saved float32 NIfTI
% beta images are compared; no double capture, no shadow wrappers, no
% global TMFC_CAP. For both estimators against the reference, over all beta files:
%   - max_abs_diff_float32, mask_mismatch_voxels_float32 (NaN pattern);
%   - max_ulp_diff_float32: exact ULP distance in float32
%     (ulp_diff_single);
%   - n_trials_nonzero and frac_trials_nonzero (trials with at least one
%     in-mask voxel differing in float32);
%   - n_voxels_nonzero, n_voxels_total and frac_voxels_nonzero (in-mask
%     voxels with a non-zero float32 difference, over all trials).
% Gates for the in-memory estimator and batched FWL solver: relative max|delta beta|/L2 <= 1e-3 and
% float32 distance <= 1 ULP, with identical masks and complete trials.
%
% Resume: a run (fir/main/in_memory/batched) with a complete output count and a
% manifest whose parent-SPM digest matches is not recomputed.
%
%   addpath('<worktree>/tests/lss_after_fir');
%   results = run_lss_after_fir_synth_tests('<out_dir>',[], '<spm_dir>', '<inputs_root>')            % all sets
%   results = run_lss_after_fir_synth_tests('<out_dir>',{{'02',1}}, '<spm_dir>', '<inputs_root>') % one set
%
% Results: <out_dir>/results.json and <out_dir>/results_table.md.
%
%
% =========================================================================
% License: GPL-3.0-or-later

if nargin < 2 || isempty(sets)
    sets = {}; ii = 0;
    for d = {'02','02v','12','15','17'}
        for k = 1:2
            ii = ii+1; sets{ii} = {d{1},k}; %#ok<AGROW>
        end
    end
end
if nargin < 4 || isempty(out_dir) || isempty(spm_dir) || isempty(inputs_root)
    error('Provide out_dir, sets, spm_dir and inputs_root.');
end

here = fileparts(mfilename('fullpath'));
addpath(spm_dir);
addpath(fullfile(here,'..','..','functions'));
addpath(fullfile(here,'legacy_main'));

main_txt = fileread(fullfile(here,'legacy_main','tmfc_LSS_after_FIR_fixed.m'));
if ~contains(main_txt,'FIR residuals are not re-scaled')
    error(['run_lss_after_fir_synth_tests: functions/tmfc_LSS_after_FIR.m is not the ' ...
           'fixed reference (residual scaling marker missing).']);
end

if ~exist(out_dir,'dir'), mkdir(out_dir); end
set(0,'DefaultFigureVisible','off');
spm('defaults','fmri');
spm_jobman('initcfg');
spm_get_defaults('cmdline',true);

results = struct();
for ic = 1:numel(sets)
    dcode = sets{ic}{1}; subj = sets{ic}{2};
    key = ['d_' dcode '_Sub_' sprintf('%02d',subj)];   % field names must not start with a digit
    fprintf('\n##### SET %s #####\n',key);
    try
        results.(key) = run_set(dcode,subj,inputs_root,out_dir);
    catch ME
        results.(key) = struct('set',sprintf('%s Sub_%02d',dcode,subj), ...
            'error',getReport(ME,'basic','hyperlinks','off'));
        fprintf('%s ERROR: %s\n',key,strtok(results.(key).error,newline));
    end
    write_case_json(results,key,out_dir);
    delete(findall(0,'type','figure'));
end

fid = fopen(fullfile(out_dir,'results.json'),'w');
fwrite(fid,jsonencode(results,'PrettyPrint',true)); fclose(fid);
write_table(results,fullfile(out_dir,'results_table.md'));
end

%% ========================================================================

function out = run_set(dcode,subj,inputs_root,out_dir)
sdir = fullfile(inputs_root,dcode,sprintf('Sub_%02d',subj));
spmmat = fullfile(sdir,'parent','SPM.mat');
if ~exist(spmmat,'file')
    error('parent SPM.mat not found: %s',spmmat);
end
sroot = fullfile(out_dir,dcode,sprintf('Sub_%02d',subj));
if ~exist(sroot,'dir'), mkdir(sroot); end
sub_name = sprintf('Sub_%02d',subj);

SPM = load(spmmat).SPM;
% conditions of interest: every condition of every session
cl = []; k = 0;
for s = 1:numel(SPM.Sess)
    for n = 1:numel(SPM.Sess(s).U)
        k = k+1;
        cl(k).sess = s; cl(k).number = n; cl(k).name = char(SPM.Sess(s).U(n).name(1));
        cl(k).file_name = ['[Run_' num2str(s) ']_[Cond_' num2str(n) ']_[' regexprep(char(SPM.Sess(s).U(n).name(1)),' ','_') ']'];
    end
end
nTrial = 0;
for k = 1:numel(cl)
    nTrial = nTrial + numel(SPM.Sess(cl(k).sess).U(cl(k).number).ons);
end

% ---- FIR task regression (once, serial) ---------------------------------
fir_root = fullfile(sroot,'fir');
mf = fullfile(sroot,'fir_manifest.json');
nRes = dir(fullfile(fir_root,'FIR_regression',sub_name,'Res_*.nii'));
if ~(exist(mf,'file') && numel(nRes) == numel(SPM.xY.VY) && manifest_ok(mf,spmmat))
    tmfc = struct();
    tmfc.project_path = fir_root;
    tmfc.subjects(1).path = spmmat; tmfc.subjects(1).name = sub_name;
    tmfc.defaults = struct('parallel',0,'maxmem',2^31,'resmem',true);
    tmfc.FIR.window = 24; tmfc.FIR.bins = 12;
    tmfc_FIR(tmfc);
    delete(findall(0,'type','figure'));
    write_manifest(mf,spmmat,struct('n_res',numel(SPM.xY.VY),'window',24,'bins',12));
end

% ---- main, in-memory estimator, batched FWL solver on the same parent and residuals ----------------------
stages = struct( ...
    'name', {'main','in_memory','batched'}, ...
    'fn',   {@tmfc_LSS_after_FIR_fixed, @tmfc_LSS_after_FIR_in_memory, ...
             @(t) tmfc_LSS_after_FIR_fwl_batch(t,1,struct('multiply','batch'))});
errs = struct(); subs = struct();
for i = 1:numel(stages)
    tag = stages(i).name;
    pdir = fullfile(sroot,tag);
    mf = fullfile(sroot,[tag '_manifest.json']);
    if exist(mf,'file') && manifest_ok(mf,spmmat) && stage_complete(pdir,sub_name,nTrial)
        errs.(tag) = '';
 subs.(tag) = [];
    else
        if exist(pdir,'dir'), rmdir(pdir,'s'); end
        mkdir(pdir);
        copyfile(fullfile(fir_root,'FIR_regression'),fullfile(pdir,'FIR_regression'));
        tmfc = struct();
        tmfc.project_path = pdir;
        tmfc.subjects(1).path = spmmat; tmfc.subjects(1).name = sub_name;
        tmfc.defaults = struct('parallel',0,'maxmem',2^31,'resmem',true);
        tmfc.session_mode = 'strict';
        tmfc.LSS_after_FIR.conditions = cl;
        err = '';
        try
            subs.(tag) = stages(i).fn(tmfc); %#ok<AGROW>
        catch ME
            err = getReport(ME,'basic','hyperlinks','off');
        end
        errs.(tag) = err; %#ok<AGROW>
        delete(findall(0,'type','figure'));
        if isempty(err)
            write_manifest(mf,spmmat,struct('n_trials',nTrial));
        end
    end
    fprintf('set %s %s: %s\n',sub_name,tag,ternary(isempty(errs.(tag)),'done','ERROR'));
end

out = struct('set',sprintf('%s Sub_%02d',dcode,subj),'nTrials',nTrial, ...
             'nScan',numel(SPM.xY.VY),'nSess',numel(SPM.Sess), ...
             'error_main',errs.main,'error_in_memory',errs.in_memory,'error_batched',errs.batched);
if ~isempty(errs.main) || ~isempty(errs.in_memory) || ~isempty(errs.batched)
    return
end

% ---- compare in-memory and batched estimators against main by the saved float32 beta files --------
mainB = fullfile(sroot,'main','LSS_regression_after_FIR',sub_name,'Betas');
out.float32.in_memory = compare_beta_dirs(mainB,fullfile(sroot,'in_memory','LSS_regression_after_FIR',sub_name,'Betas'));
out.float32.batched = compare_beta_dirs(mainB,fullfile(sroot,'batched','LSS_regression_after_FIR',sub_name,'Betas'));
out.gates.in_memory_float32_ulp_le_1 = out.float32.in_memory.max_ulp_diff_float32 <= 1 && ...
    out.float32.in_memory.mask_mismatch_voxels_float32 == 0 && out.float32.in_memory.n_compared == nTrial;
out.gates.batched_float32_ulp_le_1 = out.float32.batched.max_ulp_diff_float32 <= 1 && ...
    out.float32.batched.mask_mismatch_voxels_float32 == 0 && out.float32.batched.n_compared == nTrial;
out.gates.in_memory_rel_l2_le_1e_3 = out.float32.in_memory.max_rel_l2 <= 1e-3 && ...
    out.float32.in_memory.mask_mismatch_voxels_float32 == 0 && out.float32.in_memory.n_compared == nTrial;
out.gates.batched_rel_l2_le_1e_3 = out.float32.batched.max_rel_l2 <= 1e-3 && ...
    out.float32.batched.mask_mismatch_voxels_float32 == 0 && out.float32.batched.n_compared == nTrial;

% GLM batches, file lists, sub_check
bm = dir(fullfile(mainB,'..','GLM_batches','*.mat'));
nbe1 = 0; nbe2 = 0;
for i = 1:numel(bm)
    a = load(fullfile(bm(i).folder,bm(i).name));
    am = a.matlabbatch; am{1}.spm.stats.fmri_spec.dir = {};
    am{1}.spm.stats.fmri_spec.sess.scans = ...
        cellfun(@(p) regexprep(p,'.*[/\\]',''),am{1}.spm.stats.fmri_spec.sess.scans,'UniformOutput',false);
    c = load(fullfile(sroot,'in_memory','LSS_regression_after_FIR',sub_name,'GLM_batches',bm(i).name));
    cm = c.matlabbatch; cm{1}.spm.stats.fmri_spec.dir = {};
    cm{1}.spm.stats.fmri_spec.sess.scans = ...
        cellfun(@(p) regexprep(p,'.*[/\\]',''),cm{1}.spm.stats.fmri_spec.sess.scans,'UniformOutput',false);
    nbe1 = nbe1 + double(isequaln(am,cm));
    c = load(fullfile(sroot,'batched','LSS_regression_after_FIR',sub_name,'GLM_batches',bm(i).name));
    cm = c.matlabbatch; cm{1}.spm.stats.fmri_spec.dir = {};
    cm{1}.spm.stats.fmri_spec.sess.scans = ...
        cellfun(@(p) regexprep(p,'.*[/\\]',''),cm{1}.spm.stats.fmri_spec.sess.scans,'UniformOutput',false);
    nbe2 = nbe2 + double(isequaln(am,cm));
end
out.batches_identical_in_memory = [nbe1 numel(bm)];
out.batches_identical_batched = [nbe2 numel(bm)];
out.sub_check_identical_in_memory = ~isempty(subs.main) && ~isempty(subs.in_memory) && isequaln(subs.main,subs.in_memory);
out.sub_check_identical_batched = ~isempty(subs.main) && ~isempty(subs.batched) && isequaln(subs.main,subs.batched);
end

%% ========================================================================
%  FLOAT32 COMPARISON (saved NIfTI beta files only)
%% ========================================================================

function m = compare_beta_dirs(dirA,dirB)
% compares every Beta_*.nii of dirA with the same-named file in dirB
bb = dir(fullfile(dirA,'Beta_*.nii'));
files = {bb.name};
m = compare_beta_files(dirA,dirB,files);
end

function m = compare_beta_files(dirA,dirB,files)
% compares the named Beta_*.nii files of dirA and dirB by the saved
% float32 values: absolute difference, exact ULP distance
% (ulp_diff_single), NaN-mask mismatches and the shares of trials/voxels
% with a non-zero float32 difference
m = struct('n_compared',0,'missing_files',0,'mask_mismatch_voxels_float32',0, ...
    'max_abs_diff_float32',0,'max_ulp_diff_float32',0,'max_rel_l2',0, ...
    'n_trials_nonzero',0,'frac_trials_nonzero',NaN, ...
    'n_voxels_nonzero',0,'n_voxels_total',0,'frac_voxels_nonzero',NaN);
for i = 1:numel(files)
    fb = fullfile(dirB,files{i});
    if ~exist(fb,'file')
        m.missing_files = m.missing_files + 1;
        continue
    end
    va = spm_read_vols(spm_vol(fullfile(dirA,files{i})));
    vc = spm_read_vols(spm_vol(fb));
    m.mask_mismatch_voxels_float32 = m.mask_mismatch_voxels_float32 + ...
        nnz(isnan(va(:)) ~= isnan(vc(:)));
    b = ~isnan(va(:)) & ~isnan(vc(:));
    if ~any(b), continue; end
    a1 = va(b); c1 = vc(b);
    m.n_compared = m.n_compared + 1;
    m.n_voxels_total = m.n_voxels_total + nnz(b);
    nz = nnz(single(a1) ~= single(c1));
    m.n_voxels_nonzero = m.n_voxels_nonzero + nz;
    m.max_abs_diff_float32 = max(m.max_abs_diff_float32,max(abs(a1-c1)));
    m.max_ulp_diff_float32 = max(m.max_ulp_diff_float32,max(ulp_diff_single(a1,c1)));
    m.max_rel_l2 = max(m.max_rel_l2,max(abs(a1-c1))/max(norm(a1),eps));  % per trial: max|diff|/norm(ref)
    if nz > 0
        m.n_trials_nonzero = m.n_trials_nonzero + 1;
    end
end
if m.n_compared > 0
    m.frac_trials_nonzero = m.n_trials_nonzero / m.n_compared;
end
if m.n_voxels_total > 0
    m.frac_voxels_nonzero = m.n_voxels_nonzero / m.n_voxels_total;
end
end

function d = ulp_diff_single(a,b)
% Exact ULP distance in float32 between a and b (double or single inputs):
% single() bits via typecast, sign-magnitude -> signed int64; +0 and -0
% give 0
u1 = typecast(single(a(:)), 'uint32'); u2 = typecast(single(b(:)), 'uint32');
k1 = int64(bitand(u1, uint32(2147483647))); n1 = bitget(u1,32) ~= 0; k1(n1) = -k1(n1);
k2 = int64(bitand(u2, uint32(2147483647))); n2 = bitget(u2,32) ~= 0; k2(n2) = -k2(n2);
d = abs(k1 - k2);
end

%% ========================================================================
%  HELPERS
%% ========================================================================

function ok = manifest_ok(mf,spmmat)
ok = false;
try
    m = jsondecode(fileread(mf));
    ok = isstruct(m) && isfield(m,'sha_spmmat') && strcmp(m.sha_spmmat,sha256_file(spmmat));
catch
end
end

function write_manifest(mf,spmmat,extra)
M = struct('sha_spmmat',sha256_file(spmmat));
fn = fieldnames(extra);
for i = 1:numel(fn), M.(fn{i}) = extra.(fn{i}); end
fid = fopen(mf,'w'); fwrite(fid,jsonencode(M,'PrettyPrint',true)); fclose(fid);
end

function ok = stage_complete(pdir,sub_name,nTrial)
nb = numel(dir(fullfile(pdir,'LSS_regression_after_FIR',sub_name,'Betas','*.nii')));
ng = numel(dir(fullfile(pdir,'LSS_regression_after_FIR',sub_name,'GLM_batches','*.mat')));
ok = nb == nTrial && ng == nTrial;
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end

function h = sha256_file(fname)
% Compute a SHA-256 digest for an output file.
h = '';
try
    md  = java.security.MessageDigest.getInstance('SHA-256');
    bis = java.io.BufferedInputStream(java.io.FileInputStream(char(fname)));
    len = 65536;
    buf = zeros(1,len,'int8');
    n = bis.read(buf,0,len);
    while n > 0
        md.update(buf,0,n);
        n = bis.read(buf,0,len);
    end
    bis.close();
    d = typecast(md.digest(),'uint8');
    h = lower(sprintf('%02x',d));
catch
    [st,out] = system(['shasum -a 256 "' fname '"']);
    if st == 0
        h = strtok(out);
    else
        error('sha256_file: cannot hash %s',fname);
    end
end
end

function write_case_json(results,key,out_dir)
fid = fopen(fullfile(out_dir,['result_' key '.json']),'w');
fwrite(fid,jsonencode(results.(key),'PrettyPrint',true)); fclose(fid);
end

function write_table(results,fname)
fid = fopen(fname,'w');
fprintf(fid,'| set | nTrials | in-memory maxULP | in-memory relL2 | in-memory maxdiff f32 | in-memory mask mism | batched maxULP | batched relL2 | batched maxdiff f32 | batched mask mism | batched trials nz (frac) | batched vox nz/total (frac) | batches in-memory/batched | sub_check | gates |\n');
fprintf(fid,'|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for f = fieldnames(results)'
    r = results.(f{1});
    if ~isstruct(r) || ~isfield(r,'set'), continue; end
    if isfield(r,'error') && ~isfield(r,'nTrials')
        fprintf(fid,'| %s | error | %s |||||||||||\n',r.set,strtok(r.error,newline));
        continue
    end
    if ~isempty(r.error_main) || ~isempty(r.error_in_memory) || ~isempty(r.error_batched)
        fprintf(fid,'| %s | %d | error | main: %s | in-memory: %s | batched: %s ||||||||\n',r.set,r.nTrials, ...
            strtok(r.error_main,newline),strtok(r.error_in_memory,newline),strtok(r.error_batched,newline));
        continue
    end
    gates = 'ok';
    if ~(r.gates.in_memory_rel_l2_le_1e_3 && r.gates.batched_rel_l2_le_1e_3 && r.gates.in_memory_float32_ulp_le_1 && r.gates.batched_float32_ulp_le_1), gates = 'FAIL'; end
    fprintf(fid,'| %s | %d | %d | %.2e | %.2e | %d | %d | %.2e | %.2e | %d | %d (%.2e) | %d/%d (%.2e) | %d/%d %d/%d | %d %d | %s |\n', ...
        r.set,r.nTrials, ...
        r.float32.in_memory.max_ulp_diff_float32,r.float32.in_memory.max_rel_l2, ...
        r.float32.in_memory.max_abs_diff_float32, ...
        r.float32.in_memory.mask_mismatch_voxels_float32, ...
        r.float32.batched.max_ulp_diff_float32,r.float32.batched.max_rel_l2, ...
        r.float32.batched.max_abs_diff_float32, ...
        r.float32.batched.mask_mismatch_voxels_float32, ...
        r.float32.batched.n_trials_nonzero,r.float32.batched.frac_trials_nonzero, ...
        r.float32.batched.n_voxels_nonzero,r.float32.batched.n_voxels_total,r.float32.batched.frac_voxels_nonzero, ...
        r.batches_identical_in_memory(1),r.batches_identical_in_memory(2), ...
        r.batches_identical_batched(1),r.batches_identical_batched(2), ...
        r.sub_check_identical_in_memory,r.sub_check_identical_batched,gates);
end
fclose(fid);
end
