function [cs,fx] = peb_ppi_fixture(name,fxdir)

% Real-SPM fixtures for the tmfc_PEB_PPI tests (PEB hoist
% hoist). Builds a small estimated parent SPM (fmri_spec + fmri_est via
% spm_jobman, AR(1), canonical + time + dispersion derivatives) and the
% per-run VOI time series that tmfc_PEB_PPI reads.
%
% FORMAT [cs,fx] = peb_ppi_fixture(name,fxdir)
%   name  - 'strict' | 'repeated' | 'centering_off'
%   fxdir - disposable folder for data and the parent SPM.mat
%   cs    - case definition actually built
%   fx    - fixture bundle: spmmat, VOI series, ROI set description,
%           cond_list, session_mode, centering, expected outputs/calls
%
% Coverage:
%   strict        - session_mode 'strict', 2 runs, 2 ROIs, multi-condition
%                   sessions, bf 1/2/3 (canonical/time/dispersion), pmod
%                   column (B x Mod), 'with_mean_centering'
%   repeated      - session_mode 'repeated' (all conditions in every run),
%                   pmod, bf 1 and 2, 2 ROIs
%   centering_off - legacy cond_list without .pmod/.bf fields (backward
%                   compatibility defaults), 'no_mean_centering'
%   VOI styles    - ROI 1 "old" VOI file with field Y only; 'ROI 2' has
%                   Y + Yraw (Yraw = Y + 1000 so any accidental use of
%                   Yraw fails the numeric gate loudly); name with a space
%                   exercises regexprep in output file names
%
% STATUS: executed, PASS (see evaluation/gppi/PEB_EXACT_TESTS.md). Reproduction:
%
%   matlab -batch "cd <this_checkout>/tests/gppi/peb; results = ...
%       run_peb_ppi_tests('/tmp/tmfc10_peb_out', ...
%       '<main_toolbox>/functions', ...
%       '<spm25>'); disp(results)"
%
% License: GPL-3.0-or-later

cs = case_def(name);
fx = make_fixture(cs,fxdir);
end

%% ========================================================================

function cs = case_def(name)
cs = struct('name',name,'cvi','AR(1)','nRun',2,'seed',1, ...
            'session_mode','strict','centering','with_mean_centering', ...
            'derivs',[1 1]);
switch name
    case 'strict'
    case 'repeated',      cs.session_mode = 'repeated';      cs.seed = 2;
    case 'centering_off', cs.centering = 'no_mean_centering'; cs.seed = 3;
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

% data: task signal (incl. pmod-weighted B), drift, nuisance, AR noise
files = cell(nRun,1);  scans = cell(nRun,1); Zr = cell(nRun,1);
for r = 1:nRun
    t   = (0:nScan-1)'*TR/0.25 + 1;
    sig = zeros(nScan,4);                          % A, B, B*pmod, E
    ons = {ev{r}.ons{1} ev{r}.ons{2} ev{r}.ons{2} ev{r}.ons{3}};
    wts = {ones(numel(ev{r}.ons{1}),1) ones(numel(ev{r}.ons{2}),1) ...
          ev{r}.pmodB ones(numel(ev{r}.ons{3}),1)};  % pmod = weights only
    for c = 1:4
        u = zeros(round(nScan*TR/0.25)+80,1);
        o = round(ons{c}/0.25)+1;                   % onset times stay exact
        u(o) = wts{c};
        h = conv(u,hrf); sig(:,c) = h(t);
    end
    amp = 0.02*rand(4,nV).*(rand(4,nV)<0.4);
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
    N.mat = M; N.mat0 = M; N.descrip = 'tmfc10 peb fixture'; create(N);
    N.dat(:,:,:,:) = reshape(single(Z'),[DIM nScan]);
    files{r} = f; scans{r} = arrayfun(@(i) sprintf('%s,%d',f,i),1:nScan,'UniformOutput',false)';
end

% parent model with SPM (derivatives [1 1] so bf columns 1/2/3 exist)
pdir = fullfile(fxdir,'parent');
spec = [];
spec.dir = {fullfile(pdir)};
spec.timing = struct('units','secs','RT',TR,'fmri_t',16,'fmri_t0',8);
sess = [];
for r = 1:nRun, s = build_sess(scans{r},ev{r},nuis{r}); sess = [sess s]; end %#ok<AGROW>
spec.sess = sess;
spec.fact = struct('name',{},'levels',{});
spec.bases.hrf.derivs = cs.derivs;
spec.volt = 1;
spec.global = 'None';
spec.mthresh = 0.8;
spec.mask = {''};
spec.cvi = cs.cvi;
mb = []; mb{1}.spm.stats.fmri_spec = spec;
spm_jobman('run',mb);
spmmat = fullfile(pdir,'SPM.mat');
mb = []; mb{1}.spm.stats.fmri_est.spmmat = {spmmat};
mb{1}.spm.stats.fmri_est.write_residuals = 0;
mb{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',mb);
SPM = load(spmmat).SPM;

% VOI time series (scan time, N x 1 double): mean of small spheres
roi_names = {'ROI_1','ROI 2'};                      % second name has a space
roi_old   = [true false];                           % old VOI = Y only (no Yraw)
centers   = [5 5 4; 8 8 5];
voi = cell(nRun,numel(roi_names));
for r = 1:nRun
    for k = 1:numel(roi_names)
        d2 = (x-centers(k,1)).^2 + (y-centers(k,2)).^2 + (z-centers(k,3)).^2;
        idx = find(brain & d2 < 1.5^2);
        voi{r,k} = double(mean(Zr{r}(:,idx),2));
    end
end

% conditions of interest for tmfc_PEB_PPI
switch cs.name
    case 'strict'
        cl = [entry(SPM,1,1,1,1); entry(SPM,1,1,1,2); entry(SPM,1,2,2,1); ...
              entry(SPM,2,1,1,1); entry(SPM,2,2,1,3)];
    case 'repeated'
        cl = [entry(SPM,1,1,1,1); entry(SPM,1,2,2,1); entry(SPM,1,1,1,2)];
    case 'centering_off'                            % legacy: no .pmod/.bf fields
        cl = [legacy(SPM,1,1); legacy(SPM,1,2); legacy(SPM,2,1); legacy(SPM,2,2)];
end

% expected outputs and spm_PEB calls (mirrors tmfc_PEB_PPI loops)
nROI = numel(roi_names);
if strcmp(cs.session_mode,'repeated')
    per_run = repmat(numel(cl),1,numel(SPM.Sess));
else
    per_run = zeros(1,numel(SPM.Sess));
    for j = 1:numel(cl), per_run(cl(j).sess) = per_run(cl(j).sess) + 1; end
end
fx = struct('spmmat',spmmat,'set_name','peb_fx','roi_names',{roi_names}, ...
    'roi_old',roi_old,'voi',{voi},'conditions',cl, ...
    'session_mode',cs.session_mode,'centering',cs.centering, ...
    'expected_outputs',nROI*sum(per_run), ...
    'expected_calls_old',nROI*sum(per_run), ...
    'expected_calls_new',nROI*sum(per_run > 0));
fprintf('fixture %s: %d in-brain voxels, parent estimated: %s\n', ...
    cs.name,nnz(brain),spmmat);
end


function e = entry(SPM,s,n,p,bf)
% cond_list entry with explicit pmod column and basis function
e = struct('sess',s,'number',n,'pmod',p,'bf',bf, ...
           'name',char(SPM.Sess(s).U(n).name(p)));
end

function e = legacy(SPM,s,n)
% legacy cond_list entry (old projects): no .pmod and no .bf fields
e = struct('sess',s,'number',n,'name',char(SPM.Sess(s).U(n).name(1)));
end


function ev = make_events(T)
% Exactly 15 events: 14 random gaps rescaled to a fixed budget so the
% filter-then-index trimming of the LSS fixture pattern cannot drop any.
g = 10 + 6*rand(1,14);
g = g/sum(g) * (T - 45 - 10);                     % last onset <= T-45 < T-30
t = [10 10+cumsum(g)];
lab = [repmat([1 2],1,6) 3 3 3]; lab = lab(randperm(15));
ev.names = {'A','B','E'};
for c = 1:3
    ev.ons{c} = t(lab==c)'; ev.dur{c} = ones(nnz(lab==c),1); %#ok<AGROW>
end
ev.dur{3} = 2*ones(nnz(lab==3),1);
ev.pmodB = 0.2 + 0.8*rand(nnz(lab==2),1);          % parametric values for B
end


function s = build_sess(scans,ev,R)
s.scans = scans;
for c = 1:numel(ev.names)
    s.cond(c) = struct('name',ev.names{c},'onset',ev.ons{c},'duration',ev.dur{c}, ...
        'tmod',0,'pmod',struct('name',{},'param',{},'poly',{}),'orth',1);
end
s.cond(2).pmod = struct('name',{'Mod'},'param',{ev.pmodB},'poly',{1});
for q = 1:size(R,2)
    s.regress(q) = struct('name',sprintf('R%d',q),'val',R(:,q));
end
s.multi = {''}; s.multi_reg = {''}; s.hpf = 128;
end
