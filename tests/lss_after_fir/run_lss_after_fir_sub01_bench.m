function out = run_lss_after_fir_sub01_bench(mode,arg1,arg2,root,spmmat,spm_dir)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Benchmark the in-memory and batched LSS-after-FIR estimators against the
% fixed reference on the specified estimated parent for ds000008 sub-01:
%   <spmmat>
% All runs are serial (tmfc.defaults.parallel = 0), BLAS threads at their
% default (maxNumCompThreads is NOT called), no parfor.
%
% Comparison: ONLY the saved float32 NIfTI beta
% images are compared; there is no double capture, no shadow wrapper and
% no global TMFC_CAP, so the measured peak RSS of the runs is clean. The
% double differences are reported as metrics in the micro fixtures
% (run_lss_after_fir_micro_tests, which is not touched by this rule).
%
%   run_lss_after_fir_sub01_bench('fir',[],[],root,spmmat,spm_dir)
%       Run the FIR task regression ONCE into
%       <root>/fir/FIR_regression/sub-01/Res_*.nii (window 24 s, 12 bins).
%
%   run_lss_after_fir_sub01_bench('run',method,repeat,root,spmmat,spm_dir)
%       method = 'main' | 'in_memory' | 'batched', repeat = 1 | 2. Runs the method into
%       <root>/<method>_<repeat>/ (FIR residuals linked/copied from
%       <root>/fir) and measures the wall-clock time of the call (tic/toc).
%       The first successful 'main' run also creates the golden copy
%       <root>/golden/Betas + golden_manifest.json (SHA-256 per file).
%       in_memory/batched are compared against golden/Betas over all nTrial files
%       (Beta_[Run_r]_[Cond_c]_[name]_[Trial_n].nii from the trial maps):
%       max_abs_diff_float32, max_ulp_diff_float32 (ulp_diff_single),
%       mask_mismatch_voxels_float32, n_trials_nonzero, frac_trials_nonzero,
%       n_voxels_nonzero, n_voxels_total, frac_voxels_nonzero, n_compared.
%       Gates: in_memory gate_float32_ulp_le_1 (max_ulp_diff_float32 <= 1 &&
%       mask_mismatch_voxels_float32 == 0 && n_compared == nTrial);
%       batched gate_float32_ulp_le_1 (max_ulp_diff_float32 <= 1 &&
%       mask_mismatch_voxels_float32 == 0 && n_compared == nTrial).
%       Writes <root>/<method>_<repeat>.json.
%
%   run_lss_after_fir_sub01_bench('scaling_fix',[],[],root,spmmat,spm_dir)
%       Runs the FROZEN original routine (legacy_main/
%       tmfc_LSS_after_FIR_orig.m, which re-applies the global scaling to
%       the FIR residuals) on the same FIR residuals into <root>/orig_scaling_fix
%       (no shadow, no TMFC_CAP) and compares the saved float32 Betas with
%       golden/Betas. g and gSF per session are computed as
%       spm_fmri_spm_ui does (spm_global on the Res_*.nii; 'None':
%       gSF_r = 100/mean(g_r), 'Scaling': gSF_r = 100./g_r). Reported:
%       max_abs_diff_fixed_vs_orig_raw_float32,
%       max_abs_diff_fixed_vs_orig_over_gSF_float32,
%       max_rel_diff_fixed_vs_orig_over_gSF_float32 (voxels |va| > 1e-3
%       only, to avoid dividing by near-zero float32 noise),
%       max_ulp_diff_over_gSF_float32 (ulp_diff_single(va, single(vc/sf)))
%       and max_abs_diff_bsc_corr_float32 (beta series of session 1, first
%       condition, read from the float32 .nii files). Gates (float32):
%       identity_holds_rel_le_1e_5 and bsc_corr_invariant_le_1e_5.
%       Requires golden/Betas. Writes <root>/scaling_fix_sub01.json.
%
%   run_lss_after_fir_sub01_bench('report',[],[],root,spmmat,spm_dir)
%       Aggregates the <root>/*.json into <root>/results.json and
%       <root>/results_table.md (including the peak-RSS values recorded by
%       the bash wrapper sub01_bench_time.sh in <root>/time_logs/, when
%       present).
%
% For wall-clock AND peak RSS without cross-contamination between methods,
% run each ('run',method,repeat) in its own MATLAB process via
% sub01_bench_time.sh (/usr/bin/time -l).
%
% The LSS conditions are all conditions of all sessions of the parent
% (activation model of sub-01). Results land under
% <root>.
%
% =========================================================================
% License: GPL-3.0-or-later

if nargin < 6 || isempty(root) || isempty(spmmat) || isempty(spm_dir)
    error('Provide mode, arg1, arg2, root, spmmat and spm_dir.');
end
here    = fileparts(mfilename('fullpath'));

addpath(spm_dir);
addpath(fullfile(here,'..','..','functions'));
addpath(fullfile(here,'legacy_main'));
if ~exist(root,'dir'), mkdir(root); end
set(0,'DefaultFigureVisible','off');
spm('defaults','fmri');
spm_jobman('initcfg');
spm_get_defaults('cmdline',true);

main_txt = fileread(fullfile(here,'legacy_main','tmfc_LSS_after_FIR_fixed.m'));
if ~contains(main_txt,'FIR residuals are not re-scaled')
    error('run_lss_after_fir_sub01_bench: functions/tmfc_LSS_after_FIR.m is not the fixed reference.');
end

switch mode
    case 'fir',   out = do_fir(root,spmmat);
    case 'run',   out = do_run(root,spmmat,arg1,arg2);
    case 'scaling_fix',    out = do_scaling_fix(root,spmmat,here);
    case 'report',out = do_report(root);
    otherwise, error('unknown mode ''%s''',mode);
end
end

%% ========================================================================

function out = do_fir(root,spmmat)
fir_root = fullfile(root,'fir');
SPM = load(spmmat).SPM;
mf = fullfile(root,'fir_manifest.json');
nRes = dir(fullfile(fir_root,'FIR_regression','sub-01','Res_*.nii'));
if exist(mf,'file') && numel(nRes) == numel(SPM.xY.VY)
    out = struct('mode','fir','status','up-to-date','n_res',numel(nRes));
    return
end
tmfc = struct();
tmfc.project_path = fir_root;
tmfc.subjects(1).path = spmmat; tmfc.subjects(1).name = 'sub-01';
tmfc.defaults = struct('parallel',0,'maxmem',2^31,'resmem',true);
tmfc.FIR.window = 24; tmfc.FIR.bins = 12;
t0 = tic;
tmfc_FIR(tmfc);
delete(findall(0,'type','figure'));
nRes = dir(fullfile(fir_root,'FIR_regression','sub-01','Res_*.nii'));
write_json(mf,struct('sha_spmmat',sha256_file(spmmat),'n_res',numel(SPM.xY.VY), ...
    'wall_clock_s',toc(t0),'window',24,'bins',12,'parallel',false));
out = struct('mode','fir','status','done','n_res',numel(nRes),'wall_clock_s',toc(t0));
end

function out = do_run(root,spmmat,method,repeat)
if ~any(strcmp(method,{'main','in_memory','batched'}))
    error('method must be main, in_memory or batched');
end
if ~isnumeric(repeat) || isempty(repeat), repeat = 1; end

pdir = fullfile(root,sprintf('%s_%d',method,repeat));
if exist(pdir,'dir'), rmdir(pdir,'s'); end
mkdir(pdir);
link_or_copy(fullfile(root,'fir','FIR_regression'),fullfile(pdir,'FIR_regression'));

SPM = load(spmmat).SPM;
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

% trial maps: (run,kTrial) -> (cond,number), as in tmfc_LSS_after_FIR
maps = struct('run',{},'k',{},'cond',{},'num',{},'name',{});
for s = 1:numel(SPM.Sess)
    k = 0;
    for kc = 1:numel(cl)
        if cl(kc).sess == s
            for n = 1:numel(SPM.Sess(s).U(cl(kc).number).ons)
                k = k+1;
                maps(end+1) = struct('run',s,'k',k,'cond',cl(kc).number,'num',n, ...
                    'name',regexprep(char(SPM.Sess(s).U(cl(kc).number).name(1)),' ','_')); %#ok<AGROW>
            end
        end
    end
end
files = cell(numel(maps),1);
for t = 1:numel(maps)
    files{t} = sprintf('Beta_[Run_%d]_[Cond_%d]_[%s]_[Trial_%d].nii', ...
        maps(t).run,maps(t).cond,maps(t).name,maps(t).num);
end

tmfc = struct();
tmfc.project_path = pdir;
tmfc.subjects(1).path = spmmat; tmfc.subjects(1).name = 'sub-01';
tmfc.defaults = struct('parallel',0,'maxmem',2^31,'resmem',true);
tmfc.session_mode = 'strict';
tmfc.LSS_after_FIR.conditions = cl;

switch method
    case 'main', fn = @tmfc_LSS_after_FIR_fixed;
    case 'in_memory',   fn = @tmfc_LSS_after_FIR_in_memory;
    case 'batched',   fn = @(t) tmfc_LSS_after_FIR_fwl_batch(t,1,struct('multiply','batch'));
end

t0 = tic;
err = '';
try
    fn(tmfc);
catch ME
    err = getReport(ME,'basic','hyperlinks','off');
end
wall = toc(t0);
delete(findall(0,'type','figure'));

nb = numel(dir(fullfile(pdir,'LSS_regression_after_FIR','sub-01','Betas','*.nii')));
ng = numel(dir(fullfile(pdir,'LSS_regression_after_FIR','sub-01','GLM_batches','*.mat')));
out = struct('mode','run','method',method,'repeat',repeat,'wall_clock_s',wall, ...
    'n_trials',nTrial,'n_beta_files',nb,'n_glm_batches',ng,'error',err, ...
    'sha_spmmat',sha256_file(spmmat),'complete',nb == nTrial && ng == nTrial);

% golden: first successful main run (float32 files + manifest only)
if strcmp(method,'main') && out.complete && ~exist(fullfile(root,'golden','Betas'),'dir')
    mkdir(fullfile(root,'golden'));
    copyfile(fullfile(pdir,'LSS_regression_after_FIR','sub-01','Betas'),fullfile(root,'golden','Betas'));
    bb = dir(fullfile(root,'golden','Betas','*.nii'));
    sh = cell(numel(bb),1);
    for i = 1:numel(bb), sh{i} = struct('name',bb(i).name,'sha256',sha256_file(fullfile(bb(i).folder,bb(i).name))); end
    write_json(fullfile(root,'golden','golden_manifest.json'), ...
        struct('n_betas',numel(bb),'betas',sh,'sha_spmmat',out.sha_spmmat));
    out.is_golden = true;
end

% compare in_memory/batched against golden by the saved float32 beta files
out.float32_compare = struct('available',exist(fullfile(root,'golden','Betas'),'dir') ~= 0);
if out.float32_compare.available && out.complete && ~strcmp(method,'main')
    m = compare_beta_files(fullfile(root,'golden','Betas'), ...
        fullfile(pdir,'LSS_regression_after_FIR','sub-01','Betas'),files);
    out.float32_compare = m;
    out.float32_compare.gate_rel_l2_le_1e_3 = ...
        m.max_rel_l2 <= 1e-3 && m.mask_mismatch_voxels_float32 == 0 && m.n_compared == nTrial;
    out.float32_compare.gate_float32_ulp_le_1 = ...
        m.max_ulp_diff_float32 <= 1 && m.mask_mismatch_voxels_float32 == 0 && ...
        m.n_compared == nTrial;
end

write_json(fullfile(root,sprintf('%s_%d.json',method,repeat)),out);
fprintf('%s repeat %d: %.1f s, %d/%d betas%s\n',method,repeat,wall,nb,nTrial, ...
    ternary(isempty(err),'', [' ERROR: ' strtok(err,newline)]));
end

function out = do_scaling_fix(root,spmmat,here)
% residual scaling proof on sub-01: the FROZEN original routine, which re-applies
% gSF to the FIR residuals) vs the golden fixed main, by the saved float32
% beta files (no shadow, no TMFC_CAP)
gb = fullfile(root,'golden','Betas');
if ~isdir(gb)
    error('run_lss_after_fir_sub01_bench:scaling_fix:noGolden', ...
        'golden/Betas not found: run (''run'',''main'',1) first.');
end
addpath(fullfile(here,'legacy_main'));
cleanup = onCleanup(@() rmpath(fullfile(here,'legacy_main'))); %#ok<NASGU>

pdir = fullfile(root,'orig_scaling_fix');
if exist(pdir,'dir'), rmdir(pdir,'s'); end
mkdir(pdir);
link_or_copy(fullfile(root,'fir','FIR_regression'),fullfile(pdir,'FIR_regression'));

SPM = load(spmmat).SPM;
if size(SPM.nscan,2) ~= numel(SPM.Sess)
    error('run_lss_after_fir_sub01_bench:scaling_fix:concat','concatenated parents are not supported here.');
end
cl = []; k = 0;                                       % as in do_run
for s = 1:numel(SPM.Sess)
    for n = 1:numel(SPM.Sess(s).U)
        k = k+1;
        cl(k).sess = s; cl(k).number = n; cl(k).name = char(SPM.Sess(s).U(n).name(1));
        cl(k).file_name = ['[Run_' num2str(s) ']_[Cond_' num2str(n) ']_[' regexprep(char(SPM.Sess(s).U(n).name(1)),' ','_') ']'];
    end
end
maps = struct('run',{},'k',{},'cond',{},'num',{},'name',{});
for s = 1:numel(SPM.Sess)
    k = 0;
    for kc = 1:numel(cl)
        if cl(kc).sess == s
            for n = 1:numel(SPM.Sess(s).U(cl(kc).number).ons)
                k = k+1;
                maps(end+1) = struct('run',s,'k',k,'cond',cl(kc).number,'num',n, ...
                    'name',regexprep(char(SPM.Sess(s).U(cl(kc).number).name(1)),' ','_')); %#ok<AGROW>
            end
        end
    end
end
nTrial = numel(maps);
files = cell(nTrial,1);
for t = 1:nTrial
    files{t} = sprintf('Beta_[Run_%d]_[Cond_%d]_[%s]_[Trial_%d].nii', ...
        maps(t).run,maps(t).cond,maps(t).name,maps(t).num);
end

tmfc = struct();
tmfc.project_path = pdir;
tmfc.subjects(1).path = spmmat; tmfc.subjects(1).name = 'sub-01';
tmfc.defaults = struct('parallel',0,'maxmem',2^31,'resmem',true);
tmfc.session_mode = 'strict';
tmfc.LSS_after_FIR.conditions = cl;

err = '';
t0 = tic;
try
    tmfc_LSS_after_FIR_orig(tmfc);
catch ME
    err = getReport(ME,'basic','hyperlinks','off');
end
wall = toc(t0);
delete(findall(0,'type','figure'));
if ~isempty(err)
    out = struct('mode','scaling_fix','error',err,'wall_clock_s',wall);
    write_json(fullfile(root,'scaling_fix_sub01.json'),out);
    return
end

% g and gSF per session, as spm_fmri_spm_ui computes them on the Res_*.nii
gmode = char(SPM.xGX.iGXcalc);
isNone = strcmpi(gmode,'none');
nSess = numel(SPM.Sess);
gall = cell(1,nSess); gSF = cell(1,nSess);
for r = 1:nSess
    rows = SPM.Sess(r).row(1:SPM.nscan(r));
    gr = zeros(numel(rows),1);
    for i = 1:numel(rows)
        V = spm_data_hdr_read(fullfile(root,'fir','FIR_regression','sub-01', ...
            sprintf('Res_%04d.nii',rows(i))));
        gr(i) = spm_global(V);
    end
    gall{r} = gr;
    if isNone
        gSF{r} = 100/mean(gr)*ones(numel(rows),1);      % uniform factor per session
    else
        gSF{r} = 100./gr;
    end
end

% per-trial float32 comparison: raw, over gSF, relative (|va| > 1e-3), ULP, relL2
draw = 0; dover = 0; rel = 0; dulp = 0; rel2 = 0; nCompared = 0; mmis = 0;
for t = 1:nTrial
    if ~exist(fullfile(pdir,'LSS_regression_after_FIR','sub-01','Betas',files{t}),'file')
        continue
    end
    va = spm_read_vols(spm_vol(fullfile(gb,files{t})));
    vc = spm_read_vols(spm_vol(fullfile(pdir,'LSS_regression_after_FIR','sub-01','Betas',files{t})));
    mmis = mmis + nnz(isnan(va(:)) ~= isnan(vc(:)));
    bf = ~isnan(va(:)) & ~isnan(vc(:));
    if ~any(bf), continue; end
    a1 = va(bf); c1 = vc(bf);
    draw = max(draw,max(abs(a1-c1)));
    if isNone
        sf = gSF{maps(t).run}(1);
        over = abs(a1 - c1/sf);
        dover = max(dover,max(over));
        rel2 = max(rel2,max(over)/max(norm(a1),eps));   % per trial: max|diff|/norm(ref)
        big = abs(a1) > 1e-3;                          % skip near-zero float32 noise
        if any(big)
            rel = max(rel,max(over(big)./abs(a1(big))));
        end
        dulp = max(dulp,max(ulp_diff_single(a1,single(c1/sf))));
    end
    nCompared = nCompared + 1;
end

% BSC correlation (within session 1, first condition) from the float32 files
dr_bsc = NaN; dz_bsc = NaN; npairs = 0;
if isNone && nCompared == nTrial
    c1c = min([maps([maps.run] == 1).cond]);
    sel = find([maps.run] == 1 & [maps.cond] == c1c);
    selFiles = files(sel);
    [Bfix,Borig] = beta_series_from_files(gb, ...
        fullfile(pdir,'LSS_regression_after_FIR','sub-01','Betas'),selFiles);
    if ~isempty(Bfix) && ~isempty(Borig)
        dr_bsc = 0; dz_bsc = 0;
        nv = size(Bfix,2);
        for p = [1, 17, 101, 1001, 5001]
            if p+40 <= nv
                rf = corrcoef(Bfix(:,p),Bfix(:,p+40));
                ro = corrcoef(Borig(:,p),Borig(:,p+40));
                dr_bsc = max(dr_bsc,abs(rf(1,2)-ro(1,2)));
                dz_bsc = max(dz_bsc,dz_corr(rf(1,2),ro(1,2)));
                npairs = npairs + 1;
            end
        end
        if nv > 997
            rf = corrcoef(Bfix(:,1),Bfix(:,998));
            ro = corrcoef(Borig(:,1),Borig(:,998));
            dr_bsc = max(dr_bsc,abs(rf(1,2)-ro(1,2)));
            dz_bsc = max(dz_bsc,dz_corr(rf(1,2),ro(1,2)));
            npairs = npairs + 1;
        end
    end
end

out = struct('mode','scaling_fix','global_calculation',gmode,'nTrials',nTrial, ...
    'n_compared',[nCompared nTrial],'mask_mismatch_voxels_float32',mmis, ...
    'wall_clock_s',wall, ...
    'g_mean_per_session',[cellfun(@mean,gall)], ...
    'gSF_first_per_session',[cellfun(@(x) x(1),gSF)], ...
    'max_abs_diff_fixed_vs_orig_raw_float32',draw, ...
    'max_abs_diff_fixed_vs_orig_over_gSF_float32',dover, ...
    'max_rel_diff_fixed_vs_orig_over_gSF_float32',rel, ...
    'max_rel_l2_over_gSF_float32',rel2, ...
    'max_ulp_diff_over_gSF_float32',dulp, ...
    'max_abs_diff_bsc_corr_float32',dr_bsc, ...
    'max_abs_dz_bsc_float32',dz_bsc,'n_bsc_pairs',npairs, ...
    'identity_holds_le_1e_3',isNone && nCompared == nTrial && rel2 <= 1e-3, ...
    'bsc_dz_invariant_le_1e_3',isNone && nCompared == nTrial && dz_bsc <= 1e-3);
write_json(fullfile(root,'scaling_fix_sub01.json'),out);
fprintf('original scaling (global %s): %d/%d trials; over_gSF relL2 %.2e (ok %d), rel %.2e, ulp %d; raw %.2e; bsc dz %.2e over %d pairs (ok %d)\n', ...
    gmode,nCompared,nTrial,rel2,out.identity_holds_le_1e_3,rel,dulp,draw, ...
    dz_bsc,npairs,out.bsc_dz_invariant_le_1e_3);
end

function out = do_report(root)
jsons = dir(fullfile(root,'*.json'));
out = struct('runs',{{}});
for i = 1:numel(jsons)
    if any(strcmp(jsons(i).name,{'fir_manifest.json','results.json'})), continue; end
    try
        r = jsondecode(fileread(fullfile(jsons(i).folder,jsons(i).name)));
        if isfield(r,'method') || isfield(r,'status') || isfield(r,'mode')
            out.runs{end+1} = r; %#ok<AGROW>
        end
    catch
    end
end
% peak RSS from the bash wrapper, when present
tl = dir(fullfile(root,'time_logs','*.time'));
out.peak_rss = struct('file',{},'bytes',{});
for i = 1:numel(tl)
    txt = fileread(fullfile(tl(i).folder,tl(i).name));
    tok = regexp(txt,'(\d+)\s+maximum resident set size','tokens','once');
    if ~isempty(tok)
        out.peak_rss(end+1) = struct('file',tl(i).name,'bytes',str2double(tok{1})); %#ok<AGROW>
    end
end
% residual scaling result, when present
if exist(fullfile(root,'scaling_fix_sub01.json'),'file')
    try
        out.scaling_fix = jsondecode(fileread(fullfile(root,'scaling_fix_sub01.json')));
    catch
    end
end
write_json(fullfile(root,'results.json'),out);

fid = fopen(fullfile(root,'results_table.md'),'w');
fprintf(fid,'| run | wall (s) | complete | float32 gate | maxULP / maxdiff f32 | notes |\n');
fprintf(fid,'|---|---|---|---|---|---|\n');
for i = 1:numel(out.runs)
    r = out.runs{i};
    if isfield(r,'status') && strcmp(r.status,'up-to-date')
        fprintf(fid,'| fir | (up to date) | %d res | | | |\n',r.n_res); continue
    end
    if ~isfield(r,'method'), continue; end
    g = '-'; m = '-';
    if isfield(r,'float32_compare') && isfield(r.float32_compare,'max_ulp_diff_float32')
        m = sprintf('relL2 %.2e; %d ulp / %.2e (n %d, trials nz %d, vox nz %d/%d)', ...
            r.float32_compare.max_rel_l2, ...
            r.float32_compare.max_ulp_diff_float32,r.float32_compare.max_abs_diff_float32, ...
            r.float32_compare.n_compared,r.float32_compare.n_trials_nonzero, ...
            r.float32_compare.n_voxels_nonzero,r.float32_compare.n_voxels_total);
        if isfield(r.float32_compare,'gate_rel_l2_le_1e_3')
            g = sprintf('relL2<=1e-3 %d',r.float32_compare.gate_rel_l2_le_1e_3);
        end
        if isfield(r.float32_compare,'gate_float32_bitwise_exact')
            g = sprintf('%s; bitwise %d',g,r.float32_compare.gate_float32_bitwise_exact);
        elseif isfield(r.float32_compare,'gate_float32_ulp_le_1')
            g = sprintf('%s; ulp<=1 %d',g,r.float32_compare.gate_float32_ulp_le_1);
        end
    end
    fprintf(fid,'| %s #%d | %.1f | %d | %s | %s | %s |\n',r.method,r.repeat, ...
        r.wall_clock_s,r.complete,g,m,ternary(isempty(r.error),'',strtok(r.error,newline)));
end
for i = 1:numel(out.peak_rss)
    fprintf(fid,'| %s | peak RSS %.0f MB | | | | |\n',out.peak_rss(i).file,out.peak_rss(i).bytes/2^20);
end
if isfield(out,'scaling_fix') && isfield(out.scaling_fix,'n_compared')
    scaling_fix_result = out.scaling_fix;
    fprintf(fid,'| original scaling | %.1f | %d/%d trials, global %s | relL2 %.2e (ok %d), dz(bcs) %.2e (ok %d) | ulp %d, raw %.2e, over_gSF %.2e | %d pairs |\n', ...
        scaling_fix_result.wall_clock_s,scaling_fix_result.n_compared(1),scaling_fix_result.n_compared(2),scaling_fix_result.global_calculation, ...
        scaling_fix_result.max_rel_l2_over_gSF_float32,scaling_fix_result.identity_holds_le_1e_3, ...
        scaling_fix_result.max_abs_dz_bsc_float32,scaling_fix_result.bsc_dz_invariant_le_1e_3, ...
        scaling_fix_result.max_ulp_diff_over_gSF_float32,scaling_fix_result.max_abs_diff_fixed_vs_orig_raw_float32, ...
        scaling_fix_result.max_abs_diff_fixed_vs_orig_over_gSF_float32,scaling_fix_result.n_bsc_pairs);
end
fclose(fid);
end

%% ========================================================================
%  FLOAT32 COMPARISON HELPERS (saved NIfTI beta files only)
%% ========================================================================

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

function dz = dz_corr(r1,r2)
% |atanh(r1) - atanh(r2)| with r clamped inside (-1,1)
r1 = min(max(r1,-1+1e-12),1-1e-12);
r2 = min(max(r2,-1+1e-12),1-1e-12);
dz = abs(atanh(r1)-atanh(r2));
end

function [Vg,Vo] = beta_series_from_files(dirG,dirO,files)
% beta series (trials x voxels) read from the saved float32 beta files;
% voxels with NaN in any trial of either series are dropped. Columns are
% voxels, so corrcoef(Bfix(:,p),Bfix(:,p+off)) is the BSC correlation
% across trials between two voxels
n = numel(files);
Vg = []; Vo = [];
for i = 1:n
    vg = spm_read_vols(spm_vol(fullfile(dirG,files{i})));
    vo = spm_read_vols(spm_vol(fullfile(dirO,files{i})));
    if i == 1
        Vg = nan(numel(vg),n); Vo = nan(numel(vo),n);
    end
    Vg(:,i) = vg(:);
    Vo(:,i) = vo(:);
end
if isempty(Vg), return; end
ok = ~any(isnan(Vg),2) & ~any(isnan(Vo),2);
Vg = Vg(ok,:)'; Vo = Vo(ok,:)';
end

%% ========================================================================
%  HELPERS
%% ========================================================================

function link_or_copy(src,dst)
try
    [st] = system(sprintf('ln -s "%s" "%s"',src,dst));
    if st ~= 0, error('ln failed'); end
catch
    copyfile(src,dst);
end
end

function write_json(f,s)
fid = fopen(f,'w'); fwrite(fid,jsonencode(s,'PrettyPrint',true)); fclose(fid);
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
