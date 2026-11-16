function [cs,fx,parentinfo] = lss_rwls_fixture(name,fxdir)

% Fixtures for the LSS-rWLS tests (rWLS): synthetic runs
% with an rWLS parent, i.e. the parent SPM.mat is estimated with
% tmfc_spm_rwls_spm exactly as tmfc_rwls_LSS estimates a trial model
% (fmri_spec job with cvi AR(1), then xVi.form = 'wls' and one covariance
% component per scan over the whole parent series). Same geometry and
% machinery as lss_exact_fixture (14x14x10, TR = 2, 150 scans per run,
% 6 confound regressors, hpf 128, mthresh 0.8, AR(1) noise rho = 0.3).
%
% Cases:
%   std      - 2 runs, 15 interest trials (A/B) + 3 no-interest E events,
%              run 1 as 4D float32, run 2 as 3D int16 with per-volume scale
%   outliers - 1 run, same design, with 3 artifact scans (30, 75, 121) whose
%              noise variance is x64 in every voxel - the situation rWLS is
%              for: the per-scan variances h are strongly non-uniform
%   single   - 1 run, a single interest trial (A) + 12 no-interest E events
%              (nTrial = 1: Other_trials is the E events, non-empty)
%   concat   - 2 runs merged into one session (spm_fmri_concatenate), the
%              rWLS components span all 300 scans, K/Vi per original session
%
% FORMAT [cs,fx,parentinfo] = lss_rwls_fixture(name,fxdir)
%   name  - case name; fxdir - folder for data and the parent SPM.mat
%   cs    - case definition
%   fx    - parent SPM.mat path, LSS conditions, session mode, maxmem
%   parentinfo - form, h (per scan), hMin/hMax, nZeroW of the parent W,
%                nscan, nTrial, outlier scans
%
% License: GPL-3.0-or-later

cs = case_def(name);
[fx,parentinfo] = make_fixture(cs,fxdir);
end

%% ========================================================================

function cs = case_def(name)
cs = struct('name',name,'design','std','nRun',2,'seed',31,'outliers',[],'maxmem',2^31);
switch name
    case 'std'
    case 'outliers', cs.outliers = [30 75 121]; cs.nRun = 1; cs.seed = 32;
    case 'single',   cs.design = 'single'; cs.nRun = 1; cs.seed = 33;
    case 'concat',   cs.design = 'concat'; cs.seed = 34;
    otherwise, error('Unknown case %s',name);
end
end

function [fx,parentinfo] = make_fixture(cs,fxdir)
spm_get_defaults('stats.fmri.ufp',0.001);        % parent models with SPM defaults
if exist(fxdir,'dir'), rmdir(fxdir,'s'); end
mkdir(fxdir); mkdir(fullfile(fxdir,'data')); mkdir(fullfile(fxdir,'parent'));
rng(cs.seed);

DIM = [14 14 10]; TR = 2; nScan = 150; nRun = cs.nRun;
M   = [3 0 0 -21; 0 3 0 -21; 0 0 3 -15; 0 0 0 1];
nV  = prod(DIM);
[x,y,z] = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
r2    = ((x-7.5)/6).^2 + ((y-7.5)/6).^2 + ((z-5.5)/4.2).^2;
brain = r2 < 1;
b     = 15*rand(DIM);                                % background
b(brain) = max(800*(1+0.08*randn(nnz(brain),1)),300);
b(brain & r2 > 0.8) = 0.5*b(brain & r2 > 0.8);       % partial-volume shell (below threshold)
b = b(:)';
bi = find(brain(:) & r2(:) < 0.5); const_vox = bi(3);  % constant voxel inside the brain

% events and nuisance regressors
hrf = spm_hrf(0.25);
ev  = cell(nRun,1);  nuis = cell(nRun,1);
for r = 1:nRun
    ev{r} = make_events(cs,nScan,TR);
    nuis{r} = 0.1*cumsum(randn(nScan,6));
end

% data
files = cell(nRun,1);  scans = cell(nRun,1);
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
    drift = (0.02*randn(3,nV)).' ; basis = [cos(pi*(0:nScan-1)'/nScan) cos(2*pi*(0:nScan-1)'/nScan) (0:nScan-1)'/nScan-0.5];
    load_ = 0.002*randn(6,nV);
    Z = zeros(nScan,nV);
    for v = 1:nV
        e = filter(1,[1 -0.3],randn(nScan,1)); e = 0.015*e/std(e);
        Z(:,v) = b(v)*(1 + sig*amp(:,v) + basis*drift(v,:)' + nuis{r}*load_(:,v) + e);
    end
    % artifact scans: fresh high-variance noise in every voxel (x8 std)
    for o = cs.outliers
        Z(o,:) = Z(o,:) + b.*(0.015*8*randn(1,nV));
    end
    Z(:,const_vox) = 800;
    if r == 1, Z(10,bi(7)) = NaN; end                   % NaN voxel in one scan
    if r == 1
        f = fullfile(fxdir,'data','run1_4D.nii');          % 4D float32
        N = nifti; N.dat = file_array(f,[DIM nScan],'FLOAT32-LE',0,1,0);
        N.mat = M; N.mat0 = M; N.descrip = 'tmfc14 rwls synthetic'; create(N);
        N.dat(:,:,:,:) = reshape(single(Z'),[DIM nScan]);
        files{r} = f; scans{r} = arrayfun(@(i) sprintf('%s,%d',f,i),1:nScan,'UniformOutput',false)';
    else
        scans{r} = cell(nScan,1);                          % 3D int16 with per-volume scale
        for i = 1:nScan
            f = fullfile(fxdir,'data',sprintf('run%d_%03d.nii',r,i));
            vol = reshape(Z(i,:),DIM);
            N = nifti; N.dat = file_array(f,DIM,'INT16-LE',0,max(abs(vol(:)))/30000,0);
            N.mat = M; N.mat0 = M; create(N); N.dat(:,:,:) = vol;
            scans{r}{i} = [f ',1'];
        end
    end
end

% parent model: fmri_spec (cvi AR(1)), then the rWLS estimation
% (the pattern of tmfc_rwls_LSS, applied to the parent itself)
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
spec.global = 'None';
spec.mthresh = 0.8;
spec.mask = {''};
spec.cvi = 'AR(1)';
mb = []; mb{1}.spm.stats.fmri_spec = spec;
spm_jobman('run',mb);
spmmat = fullfile(pdir,'SPM.mat');
if strcmp(cs.design,'concat')
    spm_fmri_concatenate(spmmat,[nScan nScan]);
end
P = load(spmmat).SPM;
nScanTot = sum(P.nscan);
P.xVi.form = 'wls';
for jScan = 1:nScanTot
    P.xVi.Vi{jScan} = sparse(nScanTot,nScanTot);
    P.xVi.Vi{jScan}(jScan,jScan) = 1;
end
original_dir = pwd;
cd(pdir);
tmfc_spm_rwls_spm(P);
cd(original_dir);

% conditions of interest for LSS
SPM = load(spmmat).SPM;
names = ev{1}.names;
if strcmp(cs.design,'single'), sel = 1; else, sel = [1 2]; end
cl = []; k = 0;
for s = 1:numel(SPM.Sess)
    for n = sel
        k = k+1;
        cl(k).sess = s; cl(k).number = n; cl(k).name = names{n}; %#ok<AGROW>
        cl(k).file_name = ['[Run_' num2str(s) ']_[Cond_' num2str(n) ']_[' regexprep(char(SPM.Sess(s).U(n).name(1)),' ','_') ']']; %#ok<AGROW>
    end
end
fx = struct('maxmem',cs.maxmem,'spmmat',spmmat,'conditions',cl, ...
    'concat',strcmp(cs.design,'concat'),'nscan',SPM.nscan,'session_mode','strict');

parentinfo = struct('form',SPM.xVi.form,'h',[],'hMin',NaN,'hMax',NaN,'nZeroW',NaN, ...
    'nscan',SPM.nscan,'nTrial',sum(arrayfun(@(c) numel(SPM.Sess(cl(c).sess).U(cl(c).number).ons),1:numel(cl))), ...
    'outliers',cs.outliers);
if isfield(SPM.xVi,'h')
    parentinfo.h    = SPM.xVi.h;
    parentinfo.hMin = min(SPM.xVi.h);
    parentinfo.hMax = max(SPM.xVi.h);
end
if isfield(SPM.xX,'W')
    parentinfo.nZeroW = sum(abs(diag(full(SPM.xX.W))) <= 1e-6);
end
fprintf('fixture %s: %d in-brain voxels, rWLS parent (%d scans, h in [%g %g]): %s\n', ...
    cs.name,nnz(brain),nScanTot,parentinfo.hMin,parentinfo.hMax,spmmat);
end

function ev = make_events(cs,nScan,TR)
T = nScan*TR;
switch cs.design
    case {'std','concat'}
        t = 10 + cumsum(10+6*rand(1,16)); t = t(t < T-30);
        t = t(1:15);
        lab = [repmat([1 2],1,6) 3 3 3]; lab = lab(randperm(15));
        ev.names = {'A','B','E'};
        for c = 1:3, ev.ons{c} = t(lab==c)'; ev.dur{c} = 1*ones(nnz(lab==c),1); end
        ev.dur{3} = 2*ones(nnz(lab==3),1);
    case 'single'
        te = 20 + cumsum(6+4*rand(1,12)); te = te(te < T-30); te = te(1:12);
        ev.names = {'A','E'};
        ev.ons = {[100]',sort(te)'}; ev.dur = {[1]',2*ones(numel(te),1)}; %#ok<NBRAK1>
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

function s = build_sess(scans,ev,R)
s.scans = scans;
for c = 1:numel(ev.names)
    s.cond(c) = struct('name',ev.names{c},'onset',ev.ons{c},'duration',ev.dur{c}, ...
                       'tmod',0,'pmod',struct('name',{},'param',{},'poly',{}),'orth',1);
end
if size(R,2) > 0
    for q = 1:size(R,2), s.regress(q) = struct('name',sprintf('R%d',q),'val',R(:,q)); end
else
    s.regress = struct('name',{},'val',{});
end
s.multi = {''}; s.multi_reg = {''}; s.hpf = 128;
end
