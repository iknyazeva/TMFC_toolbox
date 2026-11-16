function results = run_lss_after_fir_micro_tests(out_dir,cases,spm_dir)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Numerical comparison of the LSS-after-FIR in-memory estimator
% (tmfc_LSS_after_FIR_in_memory) and batched FWL solver (tmfc_LSS_after_FIR_fwl_batch) with the reference
% tmfc_LSS_after_FIR on small synthetic data. One file, run from MATLAB
% (serially, BLAS threads at their default; maxNumCompThreads is NOT
% called, and neither the in-memory estimator nor the batched FWL solver uses parfor):
%
%   addpath('<worktree>/tests/lss_after_fir');
%   run_lss_after_fir_micro_tests('<out_dir>')                    % all cases
%   run_lss_after_fir_micro_tests('<out_dir>',{'std_none'})       % selected
%   run_lss_after_fir_micro_tests('<out_dir>',[],'<spm_dir>')     % SPM path
%
% For every case the script generates a fixture (synthetic data, a parent
% SPM.mat made and estimated by SPM, and the FIR task regression
% (tmfc_FIR, window 24 s, 12 bins) that creates the Res_*.nii residuals),
% then runs the fixed reference from legacy_main, the in-memory estimator,
% and the batched FWL solver on the same parent and residuals. A subset
% also exercises the trial-wise solver. The comparisons cover:
%   - trial betas in double right before the write (captured by a logging
%     wrapper of spm_data_write; one entry per chunk, concatenated to a
%     full-volume vector with NaN outside the analysis mask) and in the
%     written float32 images; gates: max|diff|/L2 <= 1e-3 and float32 distance <= 1 ULP
%     for the in-memory estimator and batched FWL solver; double max|diff| is reported as a metric;
%   - in-mask voxels (NaN pattern of the betas), NIfTI headers, the saved
%     GLM batch of every trial (isequaln), file lists and sub_check;
%   - no per-trial LSS_Run_* directories are left behind by in-memory and batched estimators.
%
% Cases: std_none, single_run, scaling, repeated, concat, single,
% single_alone (expected error in reference, in-memory and batched estimators), collinear,
% collinear_extreme (rank deficiency, batched FWL solver exact fallback), chunks (small
% stats.maxmem), const_nan (constant voxel + NaN voxel injected into the
% Res_*.nii), parallel_warning (tmfc.defaults.parallel = 1 -> warning and
% serial, bitwise equal to parallel = 0), neg_control (tampered regressor
% must be caught), plus the residual scaling proofs scaling_fix_none and scaling_fix_scaling: the fixed
% main against the original reference routine
% (tests/lss_after_fir/legacy_main/tmfc_LSS_after_FIR_orig.m) and against
% a direct OLS on the residuals:
%   scaling_fix_none:    beta_fixed = beta_orig/gSF_run (<= 1e-12) and the BSC
%               correlation (within a run) is unchanged (<= 1e-12);
%   scaling_fix_scaling: beta_fixed = direct OLS on the unscaled residuals
%               (<= 1e-10) and beta_orig = direct OLS on Y*gSF(t) (<= 1e-10).
%
% Results: <out_dir>/results.json and <out_dir>/results_table.md.
%
%
% =========================================================================
% License: GPL-3.0-or-later

if nargin < 2 || isempty(cases)
    cases = {'std_none','single_run','scaling','repeated','concat','single', ...
             'single_alone','collinear','collinear_extreme','chunks','const_nan', ...
             'parallel_warning','neg_control','scaling_fix_none','scaling_fix_scaling'};
end
if nargin < 3 || isempty(spm_dir), error('Provide out_dir, cases and spm_dir.'); end
here = fileparts(mfilename('fullpath'));
addpath(spm_dir);
addpath(fullfile(here,'..','..','functions'));
addpath(fullfile(here,'legacy_main'));

% The reference main must be the fixed reference: the
% original routine re-applies the global scaling to the FIR residuals.
main_txt = fileread(fullfile(here,'legacy_main','tmfc_LSS_after_FIR_fixed.m'));
if ~contains(main_txt,'FIR residuals are not re-scaled')
    error(['run_lss_after_fir_micro_tests: functions/tmfc_LSS_after_FIR.m is not the ' ...
           'fixed reference (residual scaling marker missing). These estimators do not re-scale ' ...
           'the FIR residuals and must be compared against the fixed main only.']);
end

if ~exist(out_dir,'dir'), mkdir(out_dir); end
set(0,'DefaultFigureVisible','off');

spm('defaults','fmri');
spm_jobman('initcfg');
spm_get_defaults('cmdline',true);

shadow_dir = fullfile(out_dir,'shadow');
make_shadow(shadow_dir,spm_dir);
addpath(shadow_dir);                       % logging wrappers first on the path
cleanup = onCleanup(@() rmpath(shadow_dir)); %#ok<NASGU>

results = struct();
for ic = 1:numel(cases)
    cs = case_def(cases{ic});
    fprintf('\n##### CASE %s #####\n',cs.name);
    switch cs.name
        case {'scaling_fix_none','scaling_fix_scaling'}
            base = 'std_none'; if strcmp(cs.name,'scaling_fix_scaling'), base = 'scaling'; end
            results = ensure_case(results,base,out_dir);
            results.(cs.name) = scaling_fix_proof(out_dir,results,cs.name);
            write_case_json(results,cs.name,out_dir);
        case 'neg_control'
            results.neg_control = neg_control(out_dir);
            write_case_json(results,'neg_control',out_dir);
        case 'parallel_warning'
            results.parallel_warning = parallel_warning(out_dir);
            write_case_json(results,'parallel_warning',out_dir);
        otherwise
            fx = make_fixture(cs,fullfile(out_dir,'fixtures',cs.name));
            [results.(cs.name),results.(['caps_' cs.name])] = run_case(cs,fx,fullfile(out_dir,'runs',cs.name));
            write_case_json(results,cs.name,out_dir);
    end
end

% the raw captured betas are not part of the report
fn = fieldnames(results);
results = rmfield(results,fn(startsWith(fn,'caps_')));

fid = fopen(fullfile(out_dir,'results.json'),'w');
fwrite(fid,jsonencode(results,'PrettyPrint',true)); fclose(fid);
write_table(results,fullfile(out_dir,'results_table.md'));
end

%% ========================================================================
%  CASE DEFINITIONS
%% ========================================================================

function cs = case_def(name)
% Parent models are always estimated with cvi AR(1) (the after-FIR LSS
% models are cvi 'None' regardless); global varies per case.
cs = struct('name',name,'global','None','mask',false, ...
            'nuis','six','design','std','nRun',2,'seed',1,'session_mode','strict', ...
            'maxmem',2^31,'inject_const_nan',false);
switch name
    case 'std_none',        % 2 runs, global 'None'
    case 'single_run',      cs.nRun = 1; cs.seed = 11;
    case 'scaling',         cs.global = 'Scaling'; cs.mask = true; cs.seed = 2;
    case 'repeated',        cs.session_mode = 'repeated'; cs.seed = 4;
    case 'concat',          cs.design = 'concat'; cs.seed = 3;
    case 'single',          cs.design = 'single';
    case 'single_alone',    cs.design = 'single_alone'; cs.nRun = 1;
    case 'collinear',       cs.design = 'collinear'; cs.nRun = 1;
    case 'collinear_extreme', cs.design = 'collinear_extreme'; cs.nRun = 1;
    case 'chunks',          cs.maxmem = 2^16; cs.seed = 5;   % ~54 voxels per chunk
    case 'const_nan',       cs.nRun = 1; cs.inject_const_nan = true; cs.seed = 6;
    case {'parallel_warning','neg_control','scaling_fix_none','scaling_fix_scaling'}
    otherwise, error('Unknown case %s',name);
end
end

%% ========================================================================
%  SHADOW WRAPPER (logging only)
%% ========================================================================

function make_shadow(shadow_dir,spm_dir)
if ~exist(shadow_dir,'dir'), mkdir(shadow_dir); end
% a copy of the SPM function under another name
txt = fileread(fullfile(spm_dir,'spm_data_write.m'));
txt = regexprep(txt,'(function[^\n]*=\s*)spm_data_write\(','$1spm_data_write_orig(','once');
fid = fopen(fullfile(shadow_dir,'spm_data_write_orig.m'),'w'); fwrite(fid,txt); fclose(fid);
w1 = {
'function V = spm_data_write(V,Y,varargin)'
'global TMFC_CAP'
'if ~isempty(TMFC_CAP) && isstruct(V) && isfield(V,''fname'') && ischar(V.fname) && ...'
'        (~isempty(strfind(V.fname,''beta_0001'')) || ~isempty(strfind(V.fname,''Beta_['')))'
'    TMFC_CAP.beta{end+1} = struct(''fname'',[V.fname ''|'' pwd],''Y'',Y(:),''nvar'',numel(varargin));'
'end'
'V = spm_data_write_orig(V,Y,varargin{:});'
};
fid = fopen(fullfile(shadow_dir,'spm_data_write.m'),'w'); fprintf(fid,'%s\n',w1{:}); fclose(fid);
end

%% ========================================================================
%  FIXTURES (parent SPM + FIR task regression)
%% ========================================================================

function fx = make_fixture(cs,fxdir)
spm_get_defaults('stats.fmri.ufp',0.001);        % parent models with SPM defaults
if exist(fxdir,'dir'), rmdir(fxdir,'s'); end
mkdir(fxdir); mkdir(fullfile(fxdir,'data')); mkdir(fullfile(fxdir,'parent')); mkdir(fullfile(fxdir,'fir'));
rng(cs.seed);

DIM = [14 14 10]; TR = 2; nScan = 150; nRun = cs.nRun;
M   = [3 0 0 -21; 0 3 0 -21; 0 0 3 -15; 0 0 0 1];
nV  = prod(DIM);
[x,y,z] = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
r2    = ((x-7.5)/6).^2 + ((y-7.5)/6).^2 + ((z-5.5)/4.2).^2;
brain = r2 < 1;
b     = 15*rand(DIM);                                % background
b(brain) = max(800*(1+0.08*randn(nnz(brain),1)),300);
b(brain & r2 > 0.8) = 0.5*b(brain & r2 > 0.8);       % partial-volume shell
b = b(:)';
bi = find(brain(:) & r2(:) < 0.5);                   % in-brain voxels

% events and nuisance regressors
hrf = spm_hrf(0.25);
ev  = cell(nRun,1);  nuis = cell(nRun,1);
for r = 1:nRun
    ev{r} = make_events(cs,r,nScan,TR);
    nuis{r} = make_nuisance(nScan);
end

% data (one 4D float32 image per run)
scans = cell(nRun,1);
for r = 1:nRun
    t   = (0:nScan-1)'*TR/0.25 + 1;
    sig = zeros(nScan,numel(ev{r}.names));
    for c = 1:numel(ev{r}.names)
        u = zeros(round(nScan*TR/0.25)+80,1);
        o = round(ev{r}.ons{c}/0.25)+1;
        u(o) = 1;
        h = conv(u,hrf); sig(:,c) = h(t);
    end
    amp = 0.02*rand(numel(ev{r}.names),nV).*(rand(numel(ev{r}.names),nV)<0.4);
    drift = (0.02*randn(3,nV)).'; basis = [cos(pi*(0:nScan-1)'/nScan) cos(2*pi*(0:nScan-1)'/nScan) (0:nScan-1)'/nScan-0.5];
    load_ = 0.002*randn(size(nuis{r},2),nV);
    Z = zeros(nScan,nV);
    for v = 1:nV
        e = filter(1,[1 -0.3],randn(nScan,1)); e = 0.015*e/std(e);
        Z(:,v) = b(v)*(1 + sig*amp(:,v) + basis*drift(v,:)' + nuis{r}*load_(:,v) + e);
    end
    Z(:,bi(3)) = 800;                                 % constant voxel inside the brain
    if r == 1, Z(10,bi(7)) = NaN; end                 % NaN voxel in one scan
    f = fullfile(fxdir,'data',sprintf('run%d_4D.nii',r));
    N = nifti; N.dat = file_array(f,[DIM nScan],'FLOAT32-LE',0,1,0);
    N.mat = M; N.mat0 = M; N.descrip = 'tmfc16 synthetic'; create(N);
    N.dat(:,:,:,:) = reshape(single(Z'),[DIM nScan]);
    scans{r} = arrayfun(@(i) sprintf('%s,%d',f,i),1:nScan,'UniformOutput',false)';
end
maskf = '';
if cs.mask
    maskf = fullfile(fxdir,'data','explicit_mask.nii');
    N = nifti; N.dat = file_array(maskf,DIM,'UINT8-LE',0,1,0); N.mat = M; N.mat0 = M; create(N);
    N.dat(:,:,:) = uint8(r2 < 0.9);
end

% parent model with SPM (cvi AR(1), as a realistic parent)
pdir = fullfile(fxdir,'parent');
spec = [];
spec.dir = {pdir};
spec.timing = struct('units','secs','RT',TR,'fmri_t',16,'fmri_t0',8);
if strcmp(cs.design,'concat')
    sess = build_sess(cat(1,scans{:}),shift_events(ev,nScan*TR),cat(1,nuis{:}));
else
    sess = [];
    for r = 1:nRun, s = build_sess(scans{r},ev{r},nuis{r}); sess = [sess s]; end %#ok<AGROW>
end
spec.sess = sess;
spec.fact = struct('name',{},'levels',{});
spec.bases.hrf.derivs = [0 0];
spec.volt = 1;
spec.global = cs.global;
spec.mthresh = 0.8;
spec.mask = {maskf};
spec.cvi = 'AR(1)';
mb = []; mb{1}.spm.stats.fmri_spec = spec;
spm_jobman('run',mb);
spmmat = fullfile(pdir,'SPM.mat');
if strcmp(cs.design,'concat')
    spm_fmri_concatenate(spmmat,[nScan nScan]);
end
mb = []; mb{1}.spm.stats.fmri_est.spmmat = {spmmat};
mb{1}.spm.stats.fmri_est.write_residuals = 0;
mb{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',mb);

% FIR task regression (creates the Res_*.nii residuals, serial)
fir_root = fullfile(fxdir,'fir');
tmfc = struct();
tmfc.project_path = fir_root;
tmfc.subjects(1).path = spmmat; tmfc.subjects(1).name = 'sub-01';
tmfc.defaults = struct('parallel',0,'maxmem',2^31,'resmem',true);
tmfc.FIR.window = 24; tmfc.FIR.bins = 12;
tmfc_FIR(tmfc);

% optional injection of a constant voxel and a NaN voxel into the residuals
inject = struct('const_vox',NaN,'nan_vox',NaN,'const_val',NaN);
if cs.inject_const_nan
    VM = spm_data_hdr_read(fullfile(pdir,'mask.nii'));
    mk = reshape(spm_data_read(VM,1:prod(DIM)),DIM) > 0;
    inmask = find(mk(:));
    inject.const_vox = inmask(2);
    inject.nan_vox   = inmask(5);
    inject.const_val = 5.0;
    rd = dir(fullfile(fir_root,'FIR_regression','sub-01','Res_*.nii'));
    for i = 1:numel(rd)
        f = fullfile(rd(i).folder,rd(i).name);
        V = spm_data_hdr_read(f);
        dat = spm_data_read(V,1:prod(DIM));
        dat(inject.const_vox) = inject.const_val;
        if i == 3, dat(inject.nan_vox) = NaN; end    % NaN in one residual volume
        V = spm_data_write(V,dat,1:prod(DIM));
    end
end

% conditions of interest for LSS
SPM = load(spmmat).SPM;
names = ev{1}.names;
switch cs.design
    case 'std',          sel = [1 2];
    case 'single',       sel = 4;
    case 'single_alone', sel = 1;
    case {'collinear','collinear_extreme'}, sel = 1;
    case 'concat',       sel = [1 2];
end
cl = []; k = 0;
nSessLSS = numel(SPM.Sess);
if any(strcmp(cs.design,{'single','single_alone','collinear','collinear_extreme'})), nSessLSS = 1; end
for s = 1:nSessLSS
    for n = sel
        k = k+1;
        cl(k).sess = s; cl(k).number = n; cl(k).name = names{n};
        cl(k).file_name = ['[Run_' num2str(s) ']_[Cond_' num2str(n) ']_[' regexprep(char(SPM.Sess(s).U(n).name(1)),' ','_') ']'];
    end
end
if strcmp(cs.session_mode,'repeated'), cl = cl([cl.sess] == 1); end
fx = struct('maxmem',cs.maxmem,'spmmat',spmmat,'conditions',cl, ...
            'concat',strcmp(cs.design,'concat'),'nscan',[nScan nScan], ...
            'session_mode',cs.session_mode,'fir_root',fir_root,'global',cs.global, ...
            'design',cs.design,'inject',inject,'dim',DIM,'n_inbrain',nnz(brain));
fprintf('fixture %s: %d in-brain voxels, parent + FIR done\n',cs.name,nnz(brain));
end

function ev = make_events(cs,r,nScan,TR)
T = nScan*TR;
switch cs.design
    case {'std','single','concat'}
        t = 10 + cumsum(10+6*rand(1,16)); t = t(t < T-30);
        t = t(1:15);
        lab = [repmat([1 2],1,6) 3 3 3]; lab = lab(randperm(15));
        ev.names = {'A','B','E'};
        for c = 1:3, ev.ons{c} = t(lab==c)'; ev.dur{c} = 1*ones(nnz(lab==c),1); end
        ev.dur{3} = 2*ones(nnz(lab==3),1);
        if strcmp(cs.design,'single') && r == 1
            ev.names{4} = 'S'; ev.ons{4} = 25.5; ev.dur{4} = 1;
        end
    case 'single_alone'
        ev.names = {'S'}; ev.ons = {100}; ev.dur = {1};
    case 'collinear'
        ev.names = {'P'}; ev.ons = {[80 80.3 200]'}; ev.dur = {[0.5 0.5 0.5]'};
    case 'collinear_extreme'
        ev.names = {'P'}; ev.ons = {[100 100.2]'}; ev.dur = {[0.5 0.5]'};
end
end

function ev = shift_events(ev,shift)
out.names = ev{1}.names;
for c = 1:numel(out.names)
    out.ons{c} = ev{1}.ons{c}; out.dur{c} = ev{1}.dur{c};
    out.ons{c} = [out.ons{c}; ev{2}.ons{c}+shift]; out.dur{c} = [out.dur{c}; ev{2}.dur{c}];
end
ev = out;
end

function R = make_nuisance(nScan)
R = 0.1*cumsum(randn(nScan,6));
end

function s = build_sess(scans,ev,R)
s.scans = scans;
for c = 1:numel(ev.names)
    s.cond(c) = struct('name',ev.names{c},'onset',ev.ons{c},'duration',ev.dur{c}, ...
                       'tmod',0,'pmod',struct('name',{},'param',{},'poly',{}),'orth',1);
end
for q = 1:size(R,2), s.regress(q) = struct('name',sprintf('R%d',q),'val',R(:,q)); end
s.multi = {''}; s.multi_reg = {''}; s.hpf = 128;
end

%% ========================================================================
%  RUN AND COMPARE
%% ========================================================================

function [out,caps] = run_case(cs,fx,run_dir)
if exist(run_dir,'dir'), rmdir(run_dir,'s'); end
mkdir(run_dir);

% run the reference, in-memory estimator and batched FWL solver (batch); batched FWL solver in trial mode on a subset
tags = {'main','in_memory','batched','integrated'};
if any(strcmp(cs.name,{'std_none','concat','collinear_extreme','chunks'}))
    tags{end+1} = 'batched_t';
end
caps = struct(); errs = struct(); eids = struct(); subs = struct(); dinfo = struct(); leak = {};
for i = 1:numel(tags)
    tag = tags{i};
    switch tag
        case 'main', fn = @tmfc_LSS_after_FIR_fixed;
        case 'in_memory',   fn = @tmfc_LSS_after_FIR_in_memory;
        case 'integrated', fn = @tmfc_LSS_after_FIR;
        case 'batched',   fn = @(t) tmfc_LSS_after_FIR_fwl_batch(t,1,struct('multiply','batch'));
        case 'batched_t',  fn = @(t) tmfc_LSS_after_FIR_fwl_batch(t,1,struct('multiply','trial'));
    end
    r = run_one(tag,fn,fx,run_dir);
    caps.(tag) = r.cap; errs.(tag) = r.err; eids.(tag) = r.errid; %#ok<AGROW>
    subs.(tag) = r.sub_check; dinfo.(tag) = r.diag; %#ok<AGROW>
    if ~strcmp(tag,'main')
        sd = dir(fullfile(run_dir,tag,'LSS_regression_after_FIR','sub-01','*'));
        leak = [leak, setdiff({sd.name},{'.','..','Betas','GLM_batches'})]; %#ok<AGROW>
    end
    delete(findall(0,'type','figure'));
end

out = struct('case',cs.name,'global',fx.global,'design',fx.design, ...
             'session_mode',fx.session_mode,'maxmem',fx.maxmem, ...
             'error_main',errs.main,'error_in_memory',errs.in_memory,'error_batched',errs.batched, ...
             'error_id_main',eids.main,'error_id_in_memory',eids.in_memory,'error_id_batched',eids.batched, ...
             'stage_extra_entries',{leak});

if strcmp(cs.name,'single_alone')
    % reference, in-memory and batched estimators must all stop with an error (empty Other_trials)
    out.single_alone_all_error = ~isempty(errs.main) && ~isempty(errs.in_memory) && ~isempty(errs.batched);
    out.in_memory_error_is_emptyOnsets = strcmp(eids.in_memory,'tmfc_LSS_after_FIR:emptyOnsets');
    out.batched_error_is_emptyOnsets = strcmp(eids.batched,'tmfc_LSS_after_FIR:emptyOnsets');
    return
end
if ~isempty(errs.main) || ~isempty(errs.in_memory) || ~isempty(errs.batched)
    return
end

% trial maps: (run,kTrial) -> (cond,number), as in tmfc_LSS_after_FIR
maps = trial_maps(fx);
nTr  = numel(maps);
out.nTrials = nTr;

% --- betas: double (captured before the write) and float32 (files) -----
keym = cellfun(@cap_key_main,caps.main.beta,'UniformOutput',false);
keystage = struct();
for i = 2:numel(tags)
    keystage.(tags{i}) = cellfun(@cap_key_stage,caps.(tags{i}).beta,'UniformOutput',false);
end
recs = struct('in_memory',[],'batched',[],'batched_t',[],'integrated',[]);
for t = 1:nTr
    km = sprintf('%d_%d',maps(t).run,maps(t).k);
    im = find(strcmp(keym,km));
    if isempty(im), continue; end
    a = cap_Y(caps.main.beta,im);
    ks = sprintf('%d_%d_%d',maps(t).run,maps(t).cond,maps(t).num);
    fn = ['Beta_[Run_' num2str(maps(t).run) ']_[Cond_' num2str(maps(t).cond) ']_[' maps(t).name ']_[Trial_' num2str(maps(t).num) '].nii'];
    va = spm_read_vols(spm_vol(fullfile(run_dir,'main','LSS_regression_after_FIR','sub-01','Betas',fn)));
    for i = 2:numel(tags)
        tag = tags{i};
        is = find(strcmp(keystage.(tag),ks));
        if isempty(is)
            recs.(tag) = [recs.(tag); struct('d',NaN,'dm',0,'f',NaN,'fm',0,'nin',0,'rel',NaN,'ulp',Inf)]; %#ok<AGROW>
            continue
        end
        c = cap_Y(caps.(tag).beta,is);
        dmask = nnz(isnan(a) ~= isnan(c));
        both = ~isnan(a) & ~isnan(c);
        dmax = max(abs(a(both)-c(both)));
        rel2 = dmax / max(norm(a(both)),eps);          % max|diff|/norm(beta_ref)
        vc = spm_read_vols(spm_vol(fullfile(run_dir,tag,'LSS_regression_after_FIR','sub-01','Betas',fn)));
        fmask = nnz(isnan(va(:)) ~= isnan(vc(:)));
        bf = ~isnan(va(:)) & ~isnan(vc(:));
        fmax = max(abs(va(bf)-vc(bf)));
        ulp = max(ulp_diff_single(va(bf),vc(bf)));
        recs.(tag) = [recs.(tag); struct('d',dmax,'dm',dmask,'f',fmax,'fm',fmask,'nin',nnz(both),'rel',rel2,'ulp',ulp)]; %#ok<AGROW>
    end
end
out.beta.in_memory = pack_beta(recs.in_memory);
out.beta.batched = pack_beta(recs.batched);
if ~isempty(recs.batched_t), out.beta.batched_t = pack_beta(recs.batched_t); end
out.gates.in_memory_float32_ulp_le_1 = out.beta.in_memory.max_ulp_float32 <= 1 && numel(recs.in_memory) == nTr;
out.gates.batched_float32_ulp_le_1 = out.beta.batched.max_ulp_float32 <= 1 && numel(recs.batched) == nTr;
out.gates.in_memory_rel_l2_le_1e_3 = out.beta.in_memory.max_rel_l2 <= 1e-3 && ...
    out.beta.in_memory.mask_mismatch_voxels_double == 0 && out.beta.in_memory.mask_mismatch_voxels_float32 == 0;
out.gates.batched_rel_l2_le_1e_3 = out.beta.batched.max_rel_l2 <= 1e-3 && ...
    out.beta.batched.mask_mismatch_voxels_double == 0 && out.beta.batched.mask_mismatch_voxels_float32 == 0;
if isfield(out.beta,'batched_t')
    out.gates.batched_t_float32_ulp_le_1 = out.beta.batched_t.max_ulp_float32 <= 1 && numel(recs.batched_t) == nTr;
    out.gates.batched_t_rel_l2_le_1e_3 = out.beta.batched_t.max_rel_l2 <= 1e-3 && ...
        out.beta.batched_t.mask_mismatch_voxels_double == 0 && out.beta.batched_t.mask_mismatch_voxels_float32 == 0;
end

% --- const_nan: the injected voxels must be NaN everywhere -------------
if cs.inject_const_nan
    fn0 = ['Beta_[Run_1]_[Cond_' num2str(maps(1).cond) ']_[' maps(1).name ']_[Trial_1].nii'];
    v0 = spm_read_vols(spm_vol(fullfile(run_dir,'main','LSS_regression_after_FIR','sub-01','Betas',fn0)));
    v1 = spm_read_vols(spm_vol(fullfile(run_dir,'in_memory','LSS_regression_after_FIR','sub-01','Betas',fn0)));
    out.const_nan.const_vox_is_nan = isnan(v0(fx.inject.const_vox)) && isnan(v1(fx.inject.const_vox));
    out.const_nan.nan_vox_is_nan   = isnan(v0(fx.inject.nan_vox)) && isnan(v1(fx.inject.nan_vox));
end

% --- batches, file lists, headers, sub_check ---------------------------
bm = dir(fullfile(run_dir,'main','LSS_regression_after_FIR','sub-01','GLM_batches','*.mat'));
out.files = struct();
nb_eq = struct('in_memory',0,'batched',0);
for i = 1:numel(bm)
    a = load(fullfile(bm(i).folder,bm(i).name));
    am = a.matlabbatch;
    am{1}.spm.stats.fmri_spec.dir = {};              % differs only by the project path
    % keep only the file name of the residual scans (the two runs have
    % different project paths)
    am{1}.spm.stats.fmri_spec.sess.scans = ...
        cellfun(@(p) regexprep(p,'.*[/\\]',''),am{1}.spm.stats.fmri_spec.sess.scans,'UniformOutput',false);
    for tag = {'in_memory','batched'}
        c = load(fullfile(run_dir,tag{1},'LSS_regression_after_FIR','sub-01','GLM_batches',bm(i).name));
        cm = c.matlabbatch;
        cm{1}.spm.stats.fmri_spec.dir = {};
        cm{1}.spm.stats.fmri_spec.sess.scans = ...
            cellfun(@(p) regexprep(p,'.*[/\\]',''),cm{1}.spm.stats.fmri_spec.sess.scans,'UniformOutput',false);
        nb_eq.(tag{1}) = nb_eq.(tag{1}) + double(isequaln(am,cm));
    end
end
out.files.batches_identical_in_memory = [nb_eq.in_memory, numel(bm)];
out.files.batches_identical_batched = [nb_eq.batched, numel(bm)];
mm = dir(fullfile(run_dir,'main','LSS_regression_after_FIR','sub-01','Betas','*.nii'));
for tag = {'in_memory','batched','batched_t'}
    ms = dir(fullfile(run_dir,tag{1},'LSS_regression_after_FIR','sub-01','Betas','*.nii'));
    out.files.(['same_beta_file_names_' tag{1}]) = isequal({mm.name},{ms.name});
end
hd = struct('in_memory',0,'batched',0);
for i = 1:numel(mm)
    a = spm_vol(fullfile(mm(i).folder,mm(i).name));
    for tag = {'in_memory','batched'}
        c = spm_vol(fullfile(run_dir,tag{1},'LSS_regression_after_FIR','sub-01','Betas',mm(i).name));
        hd.(tag{1}) = hd.(tag{1}) + ~(isequal(a.dim,c.dim) && isequal(a.mat,c.mat) && ...
            isequal(a.dt,c.dt) && isequal(a.pinfo,c.pinfo) && strcmp(a.descrip,c.descrip));
    end
end
out.files.header_mismatches = hd;
out.sub_check_identical_in_memory = isequaln(subs.main,subs.in_memory);
out.sub_check_identical_batched = isequaln(subs.main,subs.batched);

% The public entry point must match the batch estimator and preserve sub_check.
out.gates.integrated_sub_check = isequaln(subs.integrated,subs.batched);
for q = 1:numel(mm)
    vi = spm_read_vols(spm_vol(fullfile(run_dir,'integrated','LSS_regression_after_FIR','sub-01','Betas',mm(q).name)));
    vf = spm_read_vols(spm_vol(fullfile(run_dir,'batched','LSS_regression_after_FIR','sub-01','Betas',mm(q).name)));
    out.gates.integrated_sub_check = out.gates.integrated_sub_check && isequaln(vi,vf);
end
% --- diagnostics from in-memory and batched estimators --------------------------------------------
out.diag.in_memory_nMask = dinfo.in_memory(1).session(1).trial(1).nMask;
out.diag.batched_n_fallback = sum([dinfo.batched(1).session(1).trial.fallback]);
out.diag.batched_min_sv_min = min([dinfo.batched(1).session(1).trial.min_sv]);
end

function r = run_one(tag,fn,fx,run_dir)
global TMFC_CAP %#ok<GVMIS>
pdir = fullfile(run_dir,tag);
if exist(pdir,'dir'), rmdir(pdir,'s'); end
mkdir(pdir);
copyfile(fullfile(fx.fir_root,'FIR_regression'),fullfile(pdir,'FIR_regression'));
tmfc = struct();
tmfc.project_path = pdir;
tmfc.subjects(1).path = fx.spmmat; tmfc.subjects(1).name = 'sub-01';
tmfc.defaults = struct('parallel',0,'maxmem',fx.maxmem,'resmem',true);
tmfc.session_mode = fx.session_mode;
tmfc.LSS_after_FIR.conditions = fx.conditions;
TMFC_CAP = struct('beta',{{}});
err = ''; errid = ''; sc = []; di = [];
try
    if any(strcmp(tag,{'main','orig','integrated'}))
        sc = fn(tmfc);
    else
        [sc,di] = fn(tmfc);
    end
catch ME
    err = getReport(ME,'basic','hyperlinks','off');
    errid = ME.identifier;
    fprintf('%s ERROR: %s\n',tag,strtok(err,newline));
end
r = struct('cap',TMFC_CAP,'err',err,'errid',errid,'sub_check',sc,'diag',di);
TMFC_CAP = [];
end

function b = pack_beta(rec)
if isempty(rec)
    b = struct('max_abs_diff_double',NaN,'max_abs_diff_float32',NaN, ...
               'mask_mismatch_voxels_double',NaN,'mask_mismatch_voxels_float32',NaN, ...
               'n_inmask_voxels',NaN,'n_bitwise_equal_float32',NaN,'max_rel_l2',NaN,'max_ulp_float32',Inf);
    return
end
b = struct('max_abs_diff_double',max([rec.d]), ...
           'max_abs_diff_float32',max([rec.f]), ...
           'mask_mismatch_voxels_double',sum([rec.dm]), ...
           'mask_mismatch_voxels_float32',sum([rec.fm]), ...
           'n_inmask_voxels',max([rec.nin]), ...
           'n_bitwise_equal_float32',sum([rec.f] == 0), ...
           'max_rel_l2',max([rec.rel]),'max_ulp_float32',max([rec.ulp]));
end

function d = ulp_diff_single(a,b)
u1 = typecast(single(a(:)),'uint32'); u2 = typecast(single(b(:)),'uint32');
k1 = int64(bitand(u1,uint32(2147483647))); n1 = bitget(u1,32) ~= 0; k1(n1) = -k1(n1);
k2 = int64(bitand(u2,uint32(2147483647))); n2 = bitget(u2,32) ~= 0; k2(n2) = -k2(n2);
d = abs(k1-k2);
end

function k = cap_key_main(s)
if isstruct(s), s = s.fname; end
tok = regexp(s,'LSS_Run_(\d+)_Trial_(\d+)','tokens','once');
if isempty(tok), k = ''; else, k = sprintf('%s_%s',tok{1},tok{2}); end
end

function k = cap_key_stage(s)
if isstruct(s), s = s.fname; end
tok = regexp(s,'Beta_\[Run_(\d+)\]_\[Cond_(\d+)\]_\[.*\]_\[Trial_(\d+)\]','tokens','once');
if isempty(tok), k = ''; else, k = sprintf('%s_%s_%s',tok{1},tok{2},tok{3}); end
end

function v = cap_Y(betas,idx)
% concatenated captured Y of the capture entries idx (TMFC_CAP.beta is a
% cell array of structs with fields fname and Y)
e = [betas{idx}];
if isempty(e), v = []; else, v = vertcat(e.Y); end
end

%% ========================================================================
%  PARALLEL WARNING CASE
%% ========================================================================

function out = parallel_warning(out_dir)
cs = case_def('std_none');
cs.name = 'parallel_warning';
fx = make_fixture(cs,fullfile(out_dir,'fixtures','parallel_warning'));
run_dir = fullfile(out_dir,'runs','parallel_warning');
if exist(run_dir,'dir'), rmdir(run_dir,'s'); end
mkdir(run_dir);

wid = struct();
for tag = {'in_memory_p','in_memory_s','batched_p','batched_s'}
    switch tag{1}
        case {'in_memory_p','in_memory_s'}, fn = @tmfc_LSS_after_FIR_in_memory;
        otherwise,          fn = @tmfc_LSS_after_FIR;
    end
    pdir = fullfile(run_dir,tag{1});
    mkdir(pdir);
    copyfile(fullfile(fx.fir_root,'FIR_regression'),fullfile(pdir,'FIR_regression'));
    tmfc = struct();
    tmfc.project_path = pdir;
    tmfc.subjects(1).path = fx.spmmat; tmfc.subjects(1).name = 'sub-01';
    tmfc.session_mode = fx.session_mode;
    tmfc.LSS_after_FIR.conditions = fx.conditions;
    if tag{1}(end) == 'p'
        tmfc.defaults = struct('parallel',1,'maxmem',fx.maxmem,'resmem',true);
    else
        tmfc.defaults = struct('parallel',0,'maxmem',fx.maxmem,'resmem',true);
    end
    lastwarn('','');                                 % the parallel warning is the first one issued
    err = '';
    try
        fn(tmfc);
    catch ME
        err = getReport(ME,'basic','hyperlinks','off');
    end
    [wmsg,wid_tmp] = lastwarn;                       % the parallel warning is the first one issued
    wid.(tag{1}) = wid_tmp;                           %#ok<AGROW>
    wid_msg.(tag{1}) = wmsg;                          %#ok<AGROW>
    errtxt.(tag{1}) = err;                            %#ok<AGROW>
    delete(findall(0,'type','figure'));
end

% bitwise equality of the outputs: parallel=1 (warned, serial) vs parallel=0
bb = dir(fullfile(run_dir,'in_memory_s','LSS_regression_after_FIR','sub-01','Betas','*.nii'));
d1 = 0; d2 = 0;
for i = 1:numel(bb)
    va = spm_read_vols(spm_vol(fullfile(bb(i).folder,bb(i).name)));
    vb = spm_read_vols(spm_vol(fullfile(run_dir,'in_memory_p','LSS_regression_after_FIR','sub-01','Betas',bb(i).name)));
    d1 = max(d1,max(abs(va(:)-vb(:)),[],'omitnan'));
    vb = spm_read_vols(spm_vol(fullfile(run_dir,'batched_p','LSS_regression_after_FIR','sub-01','Betas',bb(i).name)));
    d2 = max(d2,max(abs(va(:)-vb(:)),[],'omitnan'));
end
out = struct('case','parallel_warning','n_beta_files',numel(bb), ...
    'warning_id_in_memory',wid.in_memory_p,'warning_id_batched',wid.batched_p, ...
    'in_memory_parallel_warning_issued',strcmp(wid.in_memory_p,'tmfc_LSS_after_FIR_in_memory:ParallelNotImplemented'), ...
    'batched_parallel_warning_issued',strcmp(wid.batched_p,'tmfc_LSS_after_FIR_fwl_batch:ParallelNotImplemented'), ...
    'in_memory_p_equals_in_memory_s_float32',d1 == 0,'batched_p_equals_in_memory_s_float32',d2 == 0, ...
    'max_abs_diff_float32_in_memory_p_vs_in_memory_s',d1,'max_abs_diff_float32_batched_p_vs_in_memory_s',d2, ...
    'error_in_memory_p',errtxt.in_memory_p,'error_batched_p',errtxt.batched_p, ...
    'warning_msg_in_memory',wid_msg.in_memory_p,'warning_msg_batched',wid_msg.batched_p);
end

%% ========================================================================
%  NEGATIVE CONTROL
%% ========================================================================

function out = neg_control(out_dir)
% Sensitivity check: shifting one onset of the first condition by 2 s in
% the parent SPM (one regressor of the design tampered) must change the
% betas well above the gates; the comparison machinery must catch it.
cs = case_def('std_none');
cs.name = 'neg_control';
fx = make_fixture(cs,fullfile(out_dir,'fixtures','neg_control'));
run_dir = fullfile(out_dir,'runs','neg_control');
if exist(run_dir,'dir'), rmdir(run_dir,'s'); end
mkdir(run_dir);

% golden: in-memory estimator and batched FWL solver on the untampered parent
run_one('in_memory',@tmfc_LSS_after_FIR_in_memory,fx,run_dir);
run_one('batched',@(t) tmfc_LSS_after_FIR_fwl_batch(t,1,struct('multiply','batch')),fx,run_dir);
delete(findall(0,'type','figure'));

% tampered parent: onset of the first A trial of run 1 shifted by 2 s
tdir = fullfile(run_dir,'tampered'); if exist(tdir,'dir'), rmdir(tdir,'s'); end
mkdir(tdir);
S = load(fx.spmmat); SPM = S.SPM;
SPM.Sess(1).U(1).ons(1) = SPM.Sess(1).U(1).ons(1) + 2;
save(fullfile(tdir,'SPM.mat'),'SPM');
fxt = fx; fxt.spmmat = fullfile(tdir,'SPM.mat');
run_one('in_memory_t',@tmfc_LSS_after_FIR_in_memory,fxt,run_dir);
run_one('batched_t',@(t) tmfc_LSS_after_FIR_fwl_batch(t,1,struct('multiply','batch')),fxt,run_dir);
delete(findall(0,'type','figure'));

bb = dir(fullfile(run_dir,'in_memory','LSS_regression_after_FIR','sub-01','Betas','*.nii'));
d1 = 0; d2 = 0;
for i = 1:numel(bb)
    va = spm_read_vols(spm_vol(fullfile(bb(i).folder,bb(i).name)));
    vb = spm_read_vols(spm_vol(fullfile(run_dir,'in_memory_t','LSS_regression_after_FIR','sub-01','Betas',bb(i).name)));
    d1 = max(d1,max(abs(va(:)-vb(:)),[],'omitnan'));
    vb = spm_read_vols(spm_vol(fullfile(run_dir,'batched_t','LSS_regression_after_FIR','sub-01','Betas',bb(i).name)));
    d2 = max(d2,max(abs(va(:)-vb(:)),[],'omitnan'));
end
out = struct('case','neg_control','negative_control',true,'n_beta_files',numel(bb), ...
    'max_abs_diff_float32_in_memory_tampered',d1,'max_abs_diff_float32_batched_tampered',d2, ...
    'tampering_detected',d1 > 1e-6 && d2 > 1e-6, ...
    'error_in_memory','','error_batched','');
end

%% ========================================================================
%  residual scaling PROOFS (fixed main vs frozen original routine vs direct OLS)
%% ========================================================================

function out = scaling_fix_proof(out_dir,results,mode)
global TMFC_CAP %#ok<GVMIS>
% runs the FROZEN original routine, which re-applies gSF to the
% residuals) on the fixture of the base case and compares it with the
% fixed main (already run in run_case) and, for 'Scaling', with a direct
% OLS on the residuals
if strcmp(mode,'scaling_fix_scaling'), base = 'scaling'; else, base = 'std_none'; end
fx = load_fixture(out_dir,base);
run_dir = fullfile(out_dir,'runs',base);
caps = results.(['caps_' base]);                     % captured double betas of main/in_memory/batched

r = run_one('orig',@(t) tmfc_LSS_after_FIR_orig(t),fx,run_dir);
delete(findall(0,'type','figure'));
if ~isempty(r.err)
    out = struct('case',mode,'error_orig',r.err);
    return
end
caps.orig = r.cap;
maps = trial_maps(fx);
nTr = numel(maps);
keym = cellfun(@cap_key_main,caps.main.beta,'UniformOutput',false);
keyo = cellfun(@cap_key_main,caps.orig.beta,'UniformOutput',false);
[~,gSF] = res_globals(fx);

if strcmp(mode,'scaling_fix_none')
    % a) beta_fixed = beta_orig / gSF_run (double, up to rounding);
    %    'None' scales the whole run by one factor, so the masks coincide
    dmax = 0; dmask = 0; bmax = 0; nin = 0; rel2 = 0;
    for t = 1:nTr
        km = sprintf('%d_%d',maps(t).run,maps(t).k);
        a = cap_Y(caps.main.beta,strcmp(keym,km));
        c = cap_Y(caps.orig.beta,strcmp(keyo,km));
        sf = gSF{maps(t).run}(1);   % 'None': one uniform factor per run
        dmask = dmask + nnz(isnan(a) ~= isnan(c));
        both = ~isnan(a) & ~isnan(c);
        dmax = max(dmax,max(abs(a(both)-c(both)/sf)));
        rel2 = max(rel2,max(abs(a(both)-c(both)/sf))/max(norm(a(both)),eps));
        bmax = max(bmax,max(abs(a(both)))); nin = max(nin,nnz(both));
    end
    % b) BSC correlation (within run 1, condition A) is unchanged by the fix
    sel = find([maps.cond] == 1 & [maps.run] == 1);
    Bfix = betaseries(caps,maps,keym,keyo,sel,'fixed');
    Borg = betaseries(caps,maps,keym,keyo,sel,'orig');
    dr = 0; dzb = 0;
    for p = [1, 5, 17]
        rf = corrcoef(Bfix(:,p),Bfix(:,p+40));
        ro = corrcoef(Borg(:,p),Borg(:,p+40));
        dr = max(dr,abs(rf(1,2)-ro(1,2)));
        dzb = max(dzb,dz_corr(rf(1,2),ro(1,2)));
    end
    out = struct('case',mode,'nTrials',nTr,'gSF_per_run',cell2mat(gSF), ...
        'max_abs_diff_beta_fixed_vs_orig_over_gSF',dmax,'max_abs_beta',bmax, ...
        'mask_mismatch_voxels',dmask,'n_inmask_voxels',nin, ...
        'max_rel_l2_over_gSF',rel2, ...
        'max_abs_diff_bsc_correlation',dr,'max_abs_dz_bsc',dzb, ...
        'identity_holds',dmax <= 1e-12,'bsc_corr_invariant',dr <= 1e-12, ...
        'gate_rel_l2_and_dz_le_1e_3',rel2 <= 1e-3 && dzb <= 1e-3);
else
    % b) direct OLS on the unscaled residuals (fixed main) and on the
    % per-scan scaled residuals Y*gSF(t) (original main); the design and
    % the data come from the job that the in-memory estimator saved in GLM_batches
    d_fix = 0; d_org = 0; mmo = 0; rel2 = 0;
    for t = 1:nTr
        a = load(fullfile(run_dir,'in_memory','LSS_regression_after_FIR','sub-01','GLM_batches', ...
            sprintf('GLM_[Run_%d]_[Cond_%d]_[%s]_[Trial_%d].mat',maps(t).run,maps(t).cond,maps(t).name,maps(t).num)));
        job = a.matlabbatch{1}.spm.stats.fmri_spec;
        if fx.concat, cc = fx.nscan; else, cc = []; end
        S = tmfc_LSS_after_FIR_design(job,cc);
        D = tmfc_LSS_after_FIR_data(job);
        XK = spm_filter(S.xX.K,S.xX.X);
        [Q,Rr] = qr(XK,0);
        km = sprintf('%d_%d',maps(t).run,maps(t).k);
        bfix = cap_Y(caps.main.beta,strcmp(keym,km));
        borg = cap_Y(caps.orig.beta,strcmp(keyo,km));
        bfix = bfix(:); borg = borg(:);                 % column vectors, as bd/bs
        bd = nan(prod(D.DIM),1); bs = nan(prod(D.DIM),1);
        sf = gSF{maps(t).run}(:);
        for i = 1:numel(D.Y)
            KY = spm_filter(S.xX.K,D.Y{i});
            tmp = Rr\(Q'*KY); bi_ = tmp(1,:);                       % direct OLS, unscaled (Current_trial row)
            tmp = Rr\(Q'*spm_filter(S.xX.K,D.Y{i}.*sf)); bs_ = tmp(1,:);  % direct OLS, scaled per scan
            bidx = D.chunk{i}(D.cmask{i});
            bd(bidx) = bi_(:);
            bs(bidx) = bs_(:);
        end
        mmo = mmo + nnz(isnan(borg) ~= isnan(bd));    % masks may differ (const voxels)
        bo = ~isnan(bfix) & ~isnan(bd);
        d_fix = max(d_fix,max(abs(bfix(bo)-bd(bo))));
        rel2 = max(rel2,max(abs(bfix(bo)-bd(bo)))/max(norm(bfix(bo)),eps));
        oo = ~isnan(borg) & ~isnan(bs);
        d_org = max(d_org,max(abs(borg(oo)-bs(oo))));
    end
    out = struct('case',mode,'nTrials',nTr, ...
        'max_abs_diff_fixed_main_vs_direct_ols_unscaled',d_fix, ...
        'max_abs_diff_orig_main_vs_direct_ols_scaled',d_org, ...
        'max_rel_l2_fixed_main_vs_direct_ols_unscaled',rel2, ...
        'mask_mismatch_orig_vs_unscaled',mmo, ...
        'fixed_equals_unscaled_ols',d_fix <= 1e-10, ...
        'orig_equals_scaled_ols',d_org <= 1e-10, ...
        'gate_rel_l2_le_1e_3',rel2 <= 1e-3);
end
TMFC_CAP = [];
end

function B = betaseries(caps,maps,keym,keyo,sel,which)
B = [];
for i = 1:numel(sel)
    t = sel(i);
    km = sprintf('%d_%d',maps(t).run,maps(t).k);
    if strcmp(which,'fixed')
        v = cap_Y(caps.main.beta,strcmp(keym,km));
    else
        v = cap_Y(caps.orig.beta,strcmp(keyo,km));
    end
    if isempty(B), B = zeros(numel(sel),numel(v)); end
    B(i,:) = v';
end
end

function dz = dz_corr(r1,r2)
% |atanh(r1) - atanh(r2)| with r clamped inside (-1,1)
r1 = min(max(r1,-1+1e-12),1-1e-12);
r2 = min(max(r2,-1+1e-12),1-1e-12);
dz = abs(atanh(r1)-atanh(r2));
end

function [g,gSF] = res_globals(fx)
% g and gSF exactly as spm_fmri_spm_ui computes them on the Res_*.nii
% ('None': one gSF per run = 100/mean(g); 'Scaling': gSF = 100./g per scan)
SPM = load(fx.spmmat).SPM;
if strcmp(fx.session_mode,'repeated'), sess_num = 1:numel(SPM.Sess); else, sess_num = unique([fx.conditions.sess]); end
g = []; gSF = cell(1,numel(sess_num));
for js = 1:numel(sess_num)
    if fx.concat
        rows = SPM.Sess(1).row(1:sum(fx.nscan));
    else
        rows = SPM.Sess(sess_num(js)).row(1:fx.nscan(min(js,numel(fx.nscan))));
    end
    gg = zeros(numel(rows),1);
    for i = 1:numel(rows)
        V = spm_data_hdr_read(fullfile(fx.fir_root,'FIR_regression','sub-01', ...
            sprintf('Res_%04d.nii',rows(i))));
        gg(i) = spm_global(V);
    end
    g = [g; gg];                                      %#ok<AGROW>
    if strcmpi(fx.global,'none')
        gSF{js} = 100/mean(gg)*ones(numel(rows),1);
    else
        gSF{js} = 100./gg;
    end
end
end

function maps = trial_maps(fx)
SPM = load(fx.spmmat).SPM;
cl = fx.conditions;
if strcmp(fx.session_mode,'repeated'), sess_num = 1:numel(SPM.Sess); else, sess_num = unique([cl.sess]); end
maps = struct('run',{},'k',{},'cond',{},'num',{},'name',{});
for js = 1:numel(sess_num)
    k = 0;
    for kc = 1:numel(cl)
        if strcmp(fx.session_mode,'repeated') || cl(kc).sess == sess_num(js)
            for n = 1:numel(SPM.Sess(sess_num(js)).U(cl(kc).number).ons)
                k = k+1;
                maps(end+1) = struct('run',sess_num(js),'k',k,'cond',cl(kc).number,'num',n, ...
                    'name',regexprep(char(SPM.Sess(sess_num(js)).U(cl(kc).number).name(1)),' ','_')); %#ok<AGROW>
            end
        end
    end
end
end

function fx = load_fixture(out_dir,name)
% fixture description rebuilt from disk (for the residual scaling proofs)
c = case_def(name);
pdir = fullfile(out_dir,'fixtures',name);
fx = struct('maxmem',c.maxmem,'spmmat',fullfile(pdir,'parent','SPM.mat'), ...
            'conditions',[],'concat',strcmp(c.design,'concat'),'nscan',[150 150], ...
            'session_mode',c.session_mode,'fir_root',fullfile(pdir,'fir'), ...
            'global',c.global,'design',c.design);
SPM = load(fx.spmmat).SPM;
sel = 1:2;
if strcmp(c.design,'single'), sel = 4; end
if any(strcmp(c.design,{'single_alone','collinear','collinear_extreme'})), sel = 1; end
cl = []; k = 0; nSessLSS = numel(SPM.Sess);
if any(strcmp(c.design,{'single','single_alone','collinear','collinear_extreme'})), nSessLSS = 1; end
for s = 1:nSessLSS
    for n = sel
        if n > numel(SPM.Sess(s).U), continue; end
        k = k+1;
        cl(k).sess = s; cl(k).number = n; cl(k).name = char(SPM.Sess(s).U(n).name(1));
        cl(k).file_name = ['[Run_' num2str(s) ']_[Cond_' num2str(n) ']_[' regexprep(char(SPM.Sess(s).U(n).name(1)),' ','_') ']'];
    end
end
if strcmp(c.session_mode,'repeated'), cl = cl([cl.sess] == 1); end
fx.conditions = cl;
end

%% ========================================================================
%  HELPERS
%% ========================================================================

function results = ensure_case(results,name,out_dir)
if ~isfield(results,name)
    cs = case_def(name);
    fx = make_fixture(cs,fullfile(out_dir,'fixtures',name));
    [results.(name),results.(['caps_' name])] = run_case(cs,fx,fullfile(out_dir,'runs',name));
    write_case_json(results,name,out_dir);
end
end

function write_case_json(results,name,out_dir)
fid = fopen(fullfile(out_dir,['result_' name '.json']),'w');
fwrite(fid,jsonencode(results.(name),'PrettyPrint',true)); fclose(fid);
end

function write_table(results,fname)
fid = fopen(fname,'w');
fprintf(fid,'| case | nTrials | in-memory double | in-memory float32 | batched double | trial-wise double | mask mism (d/f) | batches in-memory/batched | sub_check | gates |\n');
fprintf(fid,'|---|---|---|---|---|---|---|---|---|---|\n');
for f = fieldnames(results)'
    r = results.(f{1});
    if ~isstruct(r) || ~isfield(r,'case'), continue; end
    if isfield(r,'negative_control')
        fprintf(fid,'| %s | %d | tampered in-memory %.2e, batched %.2e, detected: %d ||||||||\n', ...
            r.case,r.n_beta_files,r.max_abs_diff_float32_in_memory_tampered,r.max_abs_diff_float32_batched_tampered,r.tampering_detected);
        continue
    end
    if strcmp(r.case,'parallel_warning')
        fprintf(fid,'| %s | %d | warn ids: %s / %s; in_memory_p==in_memory_s: %d, batched_p==in_memory_s: %d ||||||||\n', ...
            r.case,r.n_beta_files,r.warning_id_in_memory,r.warning_id_batched,r.in_memory_p_equals_in_memory_s_float32,r.batched_p_equals_in_memory_s_float32);
        continue
    end
    if strncmp(r.case,'scaling_fix_',numel('scaling_fix_'))
        if ~isfield(r,'nTrials')
            eo = ''; if isfield(r,'error_orig'), eo = strtok(r.error_orig,newline); end
            fprintf(fid,'| %s | error | orig: %s ||||||||\n',r.case,eo);
            continue
        end
        if isfield(r,'identity_holds')
            fprintf(fid,'| %s | %d | over_gSF %.2e, relL2 %.2e; dz(BSC) %.2e; gate(1e-3) %d; legacy 1e-12: %d/%d ||||||||\n', ...
                r.case,r.nTrials,r.max_abs_diff_beta_fixed_vs_orig_over_gSF, ...
                r.max_rel_l2_over_gSF,r.max_abs_dz_bsc,r.gate_rel_l2_and_dz_le_1e_3, ...
                r.identity_holds,r.bsc_corr_invariant);
        else
            fprintf(fid,'| %s | %d | fixed=OLS_unscaled %.2e, relL2 %.2e (gate 1e-3: %d); orig=OLS_scaled %.2e (%d); mask mism %d ||||||||\n', ...
                r.case,r.nTrials,r.max_abs_diff_fixed_main_vs_direct_ols_unscaled, ...
                r.max_rel_l2_fixed_main_vs_direct_ols_unscaled,r.gate_rel_l2_le_1e_3, ...
                r.max_abs_diff_orig_main_vs_direct_ols_scaled,r.orig_equals_scaled_ols,r.mask_mismatch_orig_vs_unscaled);
        end
        continue
    end
    if ~isempty(r.error_main) || ~isempty(r.error_in_memory) || ~isempty(r.error_batched)
        if isfield(r,'single_alone_all_error')
            fprintf(fid,'| %s | error | all error: %d; in-memory emptyOnsets: %d; batched emptyOnsets: %d ||||||||\n', ...
                r.case,r.single_alone_all_error,r.in_memory_error_is_emptyOnsets,r.batched_error_is_emptyOnsets);
        else
            fprintf(fid,'| %s | error | main: %s | in-memory: %s | batched: %s ||||||||\n', ...
                r.case,strtok(r.error_main,newline),strtok(r.error_in_memory,newline),strtok(r.error_batched,newline));
        end
        continue
    end
    batched_t = '-';
    if isfield(r.beta,'batched_t'), batched_t = sprintf('%.2e',r.beta.batched_t.max_abs_diff_double); end
    gates = 'ok';
    if ~(r.gates.in_memory_rel_l2_le_1e_3 && r.gates.batched_rel_l2_le_1e_3 && r.gates.in_memory_float32_ulp_le_1 && r.gates.batched_float32_ulp_le_1 && r.gates.integrated_sub_check), gates = 'FAIL'; end
    if isfield(r.gates,'batched_t_rel_l2_le_1e_3') && (~r.gates.batched_t_rel_l2_le_1e_3 || ~r.gates.batched_t_float32_ulp_le_1), gates = 'FAIL'; end
    fprintf(fid,'| %s | %d | %.2e | %.2e | %.2e | %s | %d / %d | %d/%d %d/%d | %d %d | %s |\n', ...
        r.case,r.nTrials,r.beta.in_memory.max_abs_diff_double,r.beta.in_memory.max_abs_diff_float32, ...
        r.beta.batched.max_abs_diff_double,batched_t, ...
        r.beta.in_memory.mask_mismatch_voxels_double,r.beta.in_memory.mask_mismatch_voxels_float32, ...
        r.files.batches_identical_in_memory(1),r.files.batches_identical_in_memory(2), ...
        r.files.batches_identical_batched(1),r.files.batches_identical_batched(2), ...
        r.sub_check_identical_in_memory,r.sub_check_identical_batched,gates);
end
fclose(fid);
end
