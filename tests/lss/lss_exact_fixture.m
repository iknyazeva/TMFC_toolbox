function [cs,fx] = lss_exact_fixture(name,fxdir)

% Synthetic fixtures for the LSS tests: case definition and parent SPM model.
% Fixed random seeds make the synthetic cases reproducible.
%
% FORMAT [cs,fx] = lss_exact_fixture(name,fxdir)
%   name  - case name; fxdir - folder for data and the parent SPM.mat
%   cs    - case definition; fx - parent SPM.mat path, LSS conditions, session mode, maxmem
%
% License: GPL-3.0-or-later

cs = case_def(name);
fx = make_fixture(cs,fxdir);
end

%% ========================================================================

function cs = case_def(name)
cs = struct('name',name,'cvi','AR(1)','global','None','mask',false, ...
            'nuis','six','design','std','nRun',2,'seed',1,'session_mode','strict','maxmem',2^31,'noconst',false);
switch name
    case 'none',             cs.cvi = 'None';
    case 'ar1'
    case 'ar1_nocv',         cs.noconst = true;          % no voxel with constant data (quasi-constant after int16 scaling)
    case 'fast',             cs.cvi = 'FAST';
    case 'ar1_scaling_mask', cs.global = 'Scaling'; cs.mask = true; cs.seed = 2;
    case 'single',           cs.design = 'single';
    case 'single_alone',     cs.design = 'single_alone'; cs.nRun = 1; cs.nuis = 'six';
    case 'collinear',        cs.design = 'collinear'; cs.nRun = 1;
    case 'collinear_extreme', cs.design = 'collinear_extreme'; cs.nRun = 1;   % two events 0.2 s apart, nothing else
    case 'rank_def',         cs.nuis = 'rankdef';
    case 'no_nuisance',      cs.nuis = 'none';
    case 'concat',           cs.design = 'concat'; cs.seed = 3;
    case 'dup_onset',        cs.design = 'dup_onset'; cs.nRun = 1;        % two events with identical onset
    case 'dense',            cs.design = 'dense'; cs.nRun = 1; cs.seed = 6;  % ISI 2.5-4.5 s, simultaneous no-interest events
    case 'repeated',         cs.session_mode = 'repeated'; cs.seed = 4;
    case 'chunks',           cs.maxmem = 2^16; cs.seed = 5;      % 54 voxels per chunk, 37 chunks
    otherwise, error('Unknown case %s',name);
end
end


function fx = make_fixture(cs,fxdir)
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
    ev{r} = make_events(cs,r,nScan,TR);
    nuis{r} = make_nuisance(cs,nScan);
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
    load_ = 0.002*randn(size(nuis{r},2),nV);
    Z = zeros(nScan,nV);
    for v = 1:nV
        e = filter(1,[1 -0.3],randn(nScan,1)); e = 0.015*e/std(e);
        Z(:,v) = b(v)*(1 + sig*amp(:,v) + basis*drift(v,:)' + nuis{r}*load_(:,v) + e);
    end
    if ~cs.noconst, Z(:,const_vox) = 800; end
    if r == 1, Z(10,bi(7)) = NaN; end                   % NaN voxel in one scan
    if r == 1
        f = fullfile(fxdir,'data','run1_4D.nii');          % 4D float32
        N = nifti; N.dat = file_array(f,[DIM nScan],'FLOAT32-LE',0,1,0);
        N.mat = M; N.mat0 = M; N.descrip = 'tmfc14 synthetic'; create(N);
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
maskf = '';
if cs.mask
    maskf = fullfile(fxdir,'data','explicit_mask.nii');
    N = nifti; N.dat = file_array(maskf,DIM,'UINT8-LE',0,1,0); N.mat = M; N.mat0 = M; create(N);
    N.dat(:,:,:) = uint8(r2 < 0.9);
end

% parent model with SPM
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
spec.cvi = cs.cvi;
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

% conditions of interest for LSS
SPM = load(spmmat).SPM;
names = ev{1}.names;
switch cs.design
    case 'std',          sel = [1 2];
    case 'single',       sel = 4;
    case 'single_alone', sel = 1;
    case {'collinear','collinear_extreme','dup_onset'}, sel = 1;
    case 'dense',        sel = [1 2];
    case 'concat',       sel = [1 2];
end
cl = []; k = 0;
nSessLSS = numel(SPM.Sess);
if any(strcmp(cs.design,{'single','single_alone','collinear','collinear_extreme','dup_onset'})), nSessLSS = 1; end
for s = 1:nSessLSS
    for n = sel
        k = k+1;
        cl(k).sess = s; cl(k).number = n; cl(k).name = names{n}; %#ok<AGROW>
        cl(k).file_name = ['[Run_' num2str(s) ']_[Cond_' num2str(n) ']_[' regexprep(char(SPM.Sess(s).U(n).name(1)),' ','_') ']']; %#ok<AGROW>
    end
end
if strcmp(cs.session_mode,'repeated'), cl = cl([cl.sess] == 1); end     % conditions of run 1 are repeated in all runs
fx = struct('maxmem',cs.maxmem,'spmmat',spmmat,'conditions',cl,'concat',strcmp(cs.design,'concat'),'nscan',[nScan nScan],'session_mode',cs.session_mode);
fprintf('fixture %s: %d in-brain voxels, parent estimated: %s\n',cs.name,nnz(brain),spmmat);
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
    case 'dup_onset'
        ev.names = {'P'}; ev.ons = {[100 100]'}; ev.dur = {[1 1]'};
    case 'dense'
        t = 10 + cumsum(2.5+2*rand(1,36)); lab = repmat([1 2],1,18); lab = lab(randperm(36));
        ev.names = {'A','B','E'};
        for c = 1:2, ev.ons{c} = t(lab==c)'; ev.dur{c} = ones(nnz(lab==c),1); end
        te = t(randperm(36,10)); ev.ons{3} = sort(te)'; ev.dur{3} = ones(10,1);   % simultaneous with interest events
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

function R = make_nuisance(cs,nScan)
switch cs.nuis
    case 'none',    R = zeros(nScan,0);
    case 'six',     R = 0.1*cumsum(randn(nScan,6));
    case 'rankdef', R = 0.1*cumsum(randn(nScan,6)); R(:,6) = R(:,1) + R(:,2); R(:,5) = 2*R(:,3);
end
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
