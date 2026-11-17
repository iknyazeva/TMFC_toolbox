function [cs,fx] = gppi_exact_fixture(name,fxdir)

% Real-SPM fixtures for the tmfc_gPPI / tmfc_gPPI_FIR exact-comparison
% tests (gPPI acceleration).
% Builds a small estimated parent SPM (fmri_spec + fmri_est via
% spm_jobman), per-run VOI time series (new Yraw+Y and old Y-only files)
% and handcrafted shared PPI regressors (zscore(psy).*zscore(VOI); NOT a
% real PEB test - PEB is validated by the separate tests/gppi/peb unit),
% plus the ROI set description.
% Pattern reused from tests/gppi/peb/peb_ppi_fixture.m (not modified).
%
% The fixture data are INPUTS: main and candidate consume identical
% files, so their outputs must match within the numeric gates of
% run_gppi_exact_tests (rel max AND rel L2 <= 1e-3, NaN masks equal,
% intentional NaN matrix diagonals). VOI Y is handcrafted: the seed
% series filtered with the parent HPF only, NOT parent-whitened and NOT
% authentic tmfc_VOI outputs - identical downstream inputs for both
% checkouts, which is all this comparison needs.
%
% Cases (first matrix = AR1 smoke across gPPI and gPPI_FIR):
%   ar1        - CORE: AR(1) parent, ROI-to-ROI (analysis 2), 2 runs,
%                3 ROIs (one old-style VOI without Yraw), multi-condition
%                sessions (A/B of interest + E nuisance), strict mode
%   fir_ar1    - same parent/VOI/PPI, FIR basis (window/bins), gPPI_FIR
%   PENDING (run only when passed explicitly; see run_gppi_exact_tests):
%   none fast rwls fir_none fir_fast fir_rwls  - remaining noise models
%   ar1_s2v    - analysis 3 (single seed-to-voxel, small analysis mask)
%   ar1_a1     - analysis 1 (seed-to-voxel + ROI-to-ROI)
%   fast_s2v fast_a1 rwls_s2v rwls_a1 and fir_ mirrors - the eight
%   FAST/rWLS x analysis 1/3 cells (completing the coverage matrix; seeds 12/13
%   per the s2v/a1 convention; executed, PASS)
%   none_s2v none_a1 fir_none_s2v fir_none_a1 - cvi 'None' via analysis
%                1/3 (policy for a1/a3: they estimate Classical via
%                spm_jobman and never call spm_est_non_sphericity; the
%                analysis-2 None failure is a documented legacy limit)
%   ar1_oldvoi - all VOIs new-style (Yraw present everywhere)
%
% STATUS: executed, PASS (evaluation/REPORT.md, Agreement). rWLS case estimates the parent with
% tmfc_spm_rwls_spm (real rWLS, xVi.form='wls'; NOT renamed AR/None).
%
% License: GPL-3.0-or-later

cs = case_def(name);
fx = make_fixture(cs,fxdir);
end

%% ========================================================================

function cs = case_def(name)
cs = struct('name',name,'cvi','AR(1)','method','gPPI','analysis',2, ...
    'old_voi',[true false false],'nRun',2,'seed',11, ...
    'window',32,'bins',8,'maxmem',2^31);
switch name
    case 'ar1'                                       % CORE
    case 'fir_ar1',   cs.method = 'fir';
    case 'none',      cs.cvi = 'None';
    case 'fast',      cs.cvi = 'FAST';
    case 'rwls',      cs.cvi = 'wls';
    case 'fir_none',  cs.method = 'fir'; cs.cvi = 'None';
    case 'fir_fast',  cs.method = 'fir'; cs.cvi = 'FAST';
    case 'fir_rwls',  cs.method = 'fir'; cs.cvi = 'wls';
    case 'ar1_s2v',   cs.analysis = 3; cs.seed = 12;
    case 'ar1_a1',    cs.analysis = 1; cs.seed = 13;
    case 'fir_ar1_s2v', cs.method = 'fir'; cs.analysis = 3; cs.seed = 12;
    case 'fir_ar1_a1',  cs.method = 'fir'; cs.analysis = 1; cs.seed = 13;
    case 'none_s2v',  cs.cvi = 'None'; cs.analysis = 3; cs.seed = 12;
    case 'none_a1',   cs.cvi = 'None'; cs.analysis = 1; cs.seed = 13;
    case 'fir_none_s2v', cs.method = 'fir'; cs.cvi = 'None'; cs.analysis = 3; cs.seed = 12;
    case 'fir_none_a1',  cs.method = 'fir'; cs.cvi = 'None'; cs.analysis = 1; cs.seed = 13;
    case 'fast_s2v',   cs.cvi = 'FAST'; cs.analysis = 3; cs.seed = 12;
    case 'fast_a1',    cs.cvi = 'FAST'; cs.analysis = 1; cs.seed = 13;
    case 'rwls_s2v',   cs.cvi = 'wls';  cs.analysis = 3; cs.seed = 12;
    case 'rwls_a1',    cs.cvi = 'wls';  cs.analysis = 1; cs.seed = 13;
    case 'fir_fast_s2v', cs.method = 'fir'; cs.cvi = 'FAST'; cs.analysis = 3; cs.seed = 12;
    case 'fir_fast_a1',  cs.method = 'fir'; cs.cvi = 'FAST'; cs.analysis = 1; cs.seed = 13;
    case 'fir_rwls_s2v', cs.method = 'fir'; cs.cvi = 'wls';  cs.analysis = 3; cs.seed = 12;
    case 'fir_rwls_a1',  cs.method = 'fir'; cs.cvi = 'wls';  cs.analysis = 1; cs.seed = 13;
    case 'ar1_oldvoi',cs.old_voi = [false false false];
    otherwise, error('Unknown case %s',name);
end
end


function fx = make_fixture(cs,fxdir)
spm_get_defaults('stats.fmri.ufp',0.001);        % parent models with SPM defaults
if exist(fxdir,'dir'), rmdir(fxdir,'s'); end
mkdir(fxdir); mkdir(fullfile(fxdir,'data')); mkdir(fullfile(fxdir,'parent'));
rng(cs.seed);

DIM = [12 12 8]; TR = 2; nScan = 120; nRun = cs.nRun;
M   = [3 0 0 -15; 0 3 0 -15; 0 0 3 -10; 0 0 0 1];
nV  = prod(DIM);
[x,y,z] = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
r2    = ((x-6.5)/5).^2 + ((y-6.5)/5).^2 + ((z-4.5)/3.5).^2;
brain = r2 < 1;
b     = 15*rand(DIM);                                % background
b(brain) = max(800*(1+0.08*randn(nnz(brain),1)),300);
b = b(:)';

% events (A/B of interest, E nuisance) and nuisance regressors
hrf  = spm_hrf(0.25);
ev   = cell(nRun,1);  nuis = cell(nRun,1);
for r = 1:nRun
    ev{r}   = make_events(nScan*TR);
    nuis{r} = 0.1*cumsum(randn(nScan,6));
end

% data: task signal, drift, nuisance, AR noise
files = cell(nRun,1);  scans = cell(nRun,1); Zr = cell(nRun,1);
for r = 1:nRun
    t   = (0:nScan-1)'*TR/0.25 + 1;
    sig = zeros(nScan,3);                           % A, B, E
    for c = 1:3
        u = zeros(round(nScan*TR/0.25)+80,1);
        o = round(ev{r}.ons{c}/0.25)+1;             % onsets stay exact
        u(o) = 1;
        h = conv(u,hrf); sig(:,c) = h(t);
    end
    amp = 0.02*rand(3,nV).*(rand(3,nV)<0.4);
    drift = (0.02*randn(3,nV)).';
    basis = [cos(pi*(0:nScan-1)'/nScan) cos(2*pi*(0:nScan-1)'/nScan) (0:nScan-1)'/nScan-0.5];
    load_ = 0.002*randn(size(nuis{r},2),nV);
    Z = zeros(nScan,nV);
    for v = 1:nV
        e = filter(1,[1 -0.3],randn(nScan,1)); e = 0.015*e/std(e);
        Z(:,v) = b(v)*(1 + sig*amp(:,v) + basis*drift(v,:)' + nuis{r}*load_(:,v) + e);
    end
    Zr{r} = Z;
    f = fullfile(fxdir,'data',sprintf('run%d_4D.nii',r));    % 4D float32
    N = nifti; N.dat = file_array(f,[DIM nScan],'FLOAT32-LE',0,1,0);
    N.mat = M; N.mat0 = M; N.descrip = 'tmfc10 gppi fixture'; create(N);
    N.dat(:,:,:,:) = reshape(single(Z'),[DIM nScan]);
    files{r} = f; scans{r} = arrayfun(@(i) sprintf('%s,%d',f,i),1:nScan,'UniformOutput',false)';
end

% parent model with SPM (spec cvi AR(1) for the rWLS parent as well)
pdir = fullfile(fxdir,'parent');
spec = [];
spec.dir = {fullfile(pdir)};
spec.timing = struct('units','secs','RT',TR,'fmri_t',16,'fmri_t0',8);
sess = [];
for r = 1:nRun, s = build_sess(scans{r},ev{r},nuis{r}); sess = [sess s]; end %#ok<AGROW>
spec.sess = sess;
spec.fact = struct('name',{},'levels',{});
spec.bases.hrf.derivs = [0 0];
spec.volt = 1;
spec.global = 'None';
spec.mthresh = 0.8;
spec.mask = {''};
% provisional valid parent ('None'); the real rWLS parent is produced
% below by tmfc_spm_rwls_spm (xVi.form='wls'); 'wls' itself is not a
% valid fmri_spec cvi menu entry. Case stays PENDING until actually run.
if strcmp(cs.cvi,'wls'), spec.cvi = 'None'; else, spec.cvi = cs.cvi; end
mb = []; mb{1}.spm.stats.fmri_spec = spec;
spm_jobman('run',mb);
spmmat = fullfile(pdir,'SPM.mat');
mb = []; mb{1}.spm.stats.fmri_est.spmmat = {spmmat};
mb{1}.spm.stats.fmri_est.write_residuals = 0;
mb{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',mb);

% rWLS parent: real wls form, one covariance component per scan
% (pattern of tests/lss/lss_rwls_fixture.m; NOT renamed AR/None).
% tmfc_spm_rwls_spm.m:350-355 raises an INTERACTIVE stop|continue
% prompt (spm_input -> questdlg) when SPM.swd already contains
% estimation files, which matlab -batch cannot show.
% Fix: re-estimate in a CLEAN sibling directory holding ONLY the
% provisional SPM.mat (design + absolute data paths): the guard never
% fires, nothing is overwritten, and the provisional estimation in
% <fxdir>/parent stays intact. The repointed spmmat below feeds the
% assert, the VOI/PPI construction and fx.spmmat (variant-neutral).
if strcmp(cs.cvi,'wls')
    pdir_wls = fullfile(fxdir,'parent_wls');
    if exist(pdir_wls,'dir'), rmdir(pdir_wls,'s'); end
    mkdir(pdir_wls);
    copyfile(spmmat,pdir_wls);          % SPM.mat only: no mask.*/beta.*/ResMS.*
    P = load(fullfile(pdir_wls,'SPM.mat')); P = P.SPM;
    P.xVi.form = 'wls';
    nScanTot = sum(P.nscan);
    for jScan = 1:nScanTot
        P.xVi.Vi{jScan} = sparse(nScanTot,nScanTot);
        P.xVi.Vi{jScan}(jScan,jScan) = 1;
    end
    P.swd = pdir_wls;                  % overwrite guard checks SPM.swd;
    original_dir = pwd;                % the helper cds to swd itself
    tmfc_spm_rwls_spm(P);
    cd(original_dir);
    spmmat = fullfile(pdir_wls,'SPM.mat');
    % the robust estimator must actually have produced a wls parent
    SPMw = load(spmmat).SPM;
    assert(strcmp(SPMw.xVi.form,'wls'), ...
        'rWLS parent xVi.form must be wls after tmfc_spm_rwls_spm');
end
SPM = load(spmmat).SPM;

% ROI set and VOI series (scan time, N x 1 double): sphere means
set_name  = 'gppi_fx'; sub = 'sub-01';
roi_names = {'ROI_1','ROI_2','ROI 3'};              % third name has a space
centers   = [5 5 4; 8 8 5; 7 4 5];
nROI      = numel(roi_names);
voi_raw   = cell(nRun,nROI); voi_y = cell(nRun,nROI);
indir = fullfile(fxdir,'inputs');
voi_dir = fullfile(indir,'ROI_sets',set_name,'VOIs',sub);
ppi_dir = fullfile(indir,'ROI_sets',set_name,'PPIs',sub);
mask_dir = fullfile(indir,'masks');
mkdir(indir); mkdir(fullfile(indir,'ROI_sets',set_name)); mkdir(voi_dir); mkdir(ppi_dir); mkdir(mask_dir);
for k = 1:nROI
    d2 = (x-centers(k,1)).^2 + (y-centers(k,2)).^2 + (z-centers(k,3)).^2;
    idx = find(brain & d2 < 1.5^2);
    N = nifti;                                       % masked ROI image
    N.dat = file_array(fullfile(mask_dir,['ROI_' regexprep(roi_names{k},' ','_') '_masked.nii']), ...
        DIM,'UINT8-LE',0,1,0);
    N.mat = M; N.mat0 = M; create(N);
    v = false(DIM); v(idx) = true; N.dat(:,:,:) = uint8(v);
    for r = 1:nRun
        voi_raw{r,k} = double(mean(Zr{r}(:,idx),2));
    end
end

% conditions of interest (strict mode, A and B in every run)
cl = [];
for r = 1:nRun
    for n = [1 2]
        nm  = char(SPM.Sess(r).U(n).name(1));
        e   = struct('sess',r,'number',n,'pmod',1,'bf',1,'name',nm, ...
            'file_name',['[Run_' num2str(r) ']_[Cond_' num2str(n) ']_[' regexprep(nm,' ','_') ']']);
        cl = [cl e]; %#ok<AGROW>
    end
end

% VOI files: Yraw + Y (new style) or Y only (old style, ROI with cs.old_voi).
% Y is handcrafted: parent-HPF-filtered seed series (not parent-whitened,
% not produced by tmfc_VOI) -- identical inputs for both checkouts
Yall = [];
for r = 1:nRun, Yall = [Yall; voi_raw{r,:}]; end %#ok<AGROW>
KYall = spm_filter(SPM.xX.K,Yall);                 % parent HPF only (handcrafted input)
for k = 1:nROI
    for r = 1:nRun
        rows = (r-1)*nScan+1:r*nScan;
        vf = fullfile(voi_dir,['VOI_' roi_names{k} '_' num2str(r) '.mat']);
        if cs.old_voi(k)
            Y = KYall(rows,k); save(vf,'Y');        % old VOI: no Yraw
        else
            Yraw = voi_raw{r,k}; Y = KYall(rows,k);
            save(vf,'Yraw','Y');
        end
    end
end

% PPI regressors: handcrafted shared inputs, zscore(psy) .* zscore(Yraw)
% per ROI x condition x run (PEB deconvolution is validated separately
% in tests/gppi/peb; identical files feed main and candidate here)
for k = 1:nROI
    for c = 1:numel(cl)
        r = cl(c).sess; n = cl(c).number;
        t = (0:nScan-1)'*TR/0.25 + 1;
        u = zeros(round(nScan*TR/0.25)+80,1);
        u(round(ev{r}.ons{n}/0.25)+1) = 1;
        psy = conv(u,hrf); psy = psy(t);
        PPI.name = ['PPI_' cl(c).name roi_names{k}];
        PPI.ppi  = zscore(psy).*zscore(voi_raw{r,k});
        pf = fullfile(ppi_dir,['PPI_[' regexprep(roi_names{k},' ','_') ']_' cl(c).file_name '.mat']);
        save(pf,'PPI');
    end
end

if strcmp(cs.method,'fir'), fxfunc = 'tmfc_gPPI_FIR'; meth_dir = 'gPPI_FIR';
else, fxfunc = 'tmfc_gPPI'; meth_dir = 'gPPI'; end
fx = struct('spmmat',spmmat,'set_name',set_name,'roi_names',{roi_names}, ...
    'old_voi',cs.old_voi,'conditions',cl,'session_mode','strict', ...
    'method',cs.method,'func',fxfunc,'meth_dir',meth_dir,'analysis',cs.analysis, ...
    'window',cs.window,'bins',cs.bins,'maxmem',cs.maxmem,'inputs_dir',indir, ...
    'subject',sub,'nRun',nRun,'cvi',cs.cvi, ...
    'expected_batches',nROI,'expected_matrices',numel(cl), ...
    'expected_betas_per_roi',numel(cl));
fprintf('fixture %s: %d in-brain voxels, parent estimated: %s\n', ...
    cs.name,nnz(brain),spmmat);
end


function ev = make_events(T)
% Exactly 15 events with rescaled gaps (LSS/PEB fixture pattern)
g = 10 + 6*rand(1,14);
g = g/sum(g) * (T - 45 - 10);
t = [10 10+cumsum(g)];
lab = [repmat([1 2],1,6) 3 3 3]; lab = lab(randperm(15));
ev.names = {'A','B','E'};
for c = 1:3
    ev.ons{c} = t(lab==c)'; ev.dur{c} = ones(nnz(lab==c),1); %#ok<AGROW>
end
ev.dur{3} = 2*ones(nnz(lab==3),1);
end


function s = build_sess(scans,ev,R)
s.scans = scans;
for c = 1:numel(ev.names)
    s.cond(c) = struct('name',ev.names{c},'onset',ev.ons{c},'duration',ev.dur{c}, ...
        'tmod',0,'pmod',struct('name',{},'param',{},'poly',{}),'orth',1);
end
for q = 1:size(R,2)
    s.regress(q) = struct('name',sprintf('R%d',q),'val',R(:,q));
end
s.multi = {''}; s.multi_reg = {''}; s.hpf = 128;
end
