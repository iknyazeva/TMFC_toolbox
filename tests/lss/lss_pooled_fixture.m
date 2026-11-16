function [cs,fx] = lss_pooled_fixture(name,fxdir)

% Synthetic fixtures for the LSS pooled tests: runs in which the per-trial
% whitening W_k genuinely differs between trials. Same geometry and machinery
% as lss_exact_fixture (14x14x10, TR = 2, 150 scans, 6 confound regressors,
% AR(1) parent, hpf 128, mthresh 0.8), with
%
%   - AR(1) noise with rho = 0.5 (instead of 0.3);
%   - condition-dependent amplitudes (A strong 0.06, B weak 0.015) and a
%     log-normal per-trial amplitude spread 10^1.2*(2u-1) (~ x0.06 .. x16),
%     so the lumped Other_trials column fits the remaining trials badly and
%     ResSS_k (and therefore Cy_k and W_k) varies strongly over trials;
%   - E as a no-interest condition (leakage into the Other_trials residuals).
%
% Cases:
%   var_std   - sparse design: 15 interest trials, ISI 10-16 s, duration 1 s
%   var_dense - dense design: 36 interest trials, ISI 2.5-4.5 s, duration 1 s,
%               plus 10 no-interest E events sharing onsets with interest trials
%
% FORMAT [cs,fx] = lss_pooled_fixture(name,fxdir)
%   name - 'var_std' | 'var_dense'; fxdir - folder for data and the parent SPM.mat
%   cs   - case definition; fx - parent SPM.mat path, LSS conditions, session mode
%
% License: GPL-3.0-or-later

switch name
    case 'var_std',   seed = 21; isi = [10 6]; nEv = 15;
    case 'var_dense', seed = 22; isi = [2.5 2]; nEv = 36;
    otherwise, error('Unknown case %s',name);
end
cs = struct('name',name,'cvi','AR(1)','seed',seed,'maxmem',2^31);

spm_get_defaults('stats.fmri.ufp',0.001);        % parent models with SPM defaults
if exist(fxdir,'dir'), rmdir(fxdir,'s'); end
mkdir(fxdir); mkdir(fullfile(fxdir,'data')); mkdir(fullfile(fxdir,'parent'));
rng(seed);

DIM = [14 14 10]; TR = 2; nScan = 150;
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

% events: interest trials of conditions A (strong) and B (weak), E = no interest
t   = 10 + cumsum(isi(1)+isi(2)*rand(1,nEv+5)); t = t(t < nScan*TR-30);
t   = t(1:nEv);
lab = repmat([1 2],1,ceil(nEv/2)); lab = lab(1:nEv); lab = lab(randperm(nEv));
onsA = t(lab==1)'; onsB = t(lab==2)';
te   = t(randperm(nEv,10)); onsE = sort(te)';        % E shares onsets with interest trials
amp0 = [0.06 0.015 0.01];                            % A, B, E
mlt  = 10.^(1.2*(2*rand(nEv,1)-1));                  % per-trial lognormal spread

% nuisance regressors and drift
nuis  = 0.1*cumsum(randn(nScan,6));
drift = (0.02*randn(3,nV)).';
basis = [cos(pi*(0:nScan-1)'/nScan) cos(2*pi*(0:nScan-1)'/nScan) (0:nScan-1)'/nScan-0.5];
load_ = 0.002*randn(6,nV);

% convolved regressors: one column per interest trial + one for E
hrf  = spm_hrf(0.25);
res_ = 0.25; len = round(nScan*TR/res_);
tt   = (0:nScan-1)'*TR/res_ + 1;
onsv = [onsA; onsB]; durv = [ones(numel(onsA)+numel(onsB),1); 2*ones(numel(onsE),1)];
condOfTrial = [ones(numel(onsA),1); 2*ones(numel(onsB),1)];
nTrial = numel(onsv);
sig = zeros(nScan,nTrial+1);
for c = 1:nTrial+1
    u = zeros(len+80,1);
    if c <= nTrial, o = onsv(c); d = 1; else, o = onsE; d = durv(end); end
    for k = 1:numel(o)
        a = round(o(k)/res_)+1; bb = round((o(k)+d)/res_);
        if a <= len, u(a:min(max(bb,a),len)) = 1; end
    end
    h = conv(u,hrf); sig(:,c) = h(tt);
end

% spatial amplitudes: condition pattern x per-trial multiplier
base = 0.35*rand(2,nV);
amp  = zeros(nTrial+1,nV);
for k = 1:nTrial
    amp(k,:) = amp0(condOfTrial(k))*mlt(k)*base(condOfTrial(k),:);
end
amp(end,:) = amp0(3)*0.5*rand(1,nV);

% data
Z = zeros(nScan,nV);
for v = 1:nV
    e = filter(1,[1 -0.5],randn(nScan,1)); e = 0.015*e/std(e);   % AR(1), rho = 0.5
    Z(:,v) = b(v)*(1 + sig*amp(:,v) + basis*drift(v,:)' + nuis*load_(:,v) + e);
end
Z(:,const_vox) = 800;
Z(10,bi(7)) = NaN;                                    % NaN voxel in one scan

f = fullfile(fxdir,'data','run1_4D.nii');             % 4D float32
N = nifti; N.dat = file_array(f,[DIM nScan],'FLOAT32-LE',0,1,0);
N.mat = M; N.mat0 = M; N.descrip = 'tmfc14 s3 synthetic'; create(N);
N.dat(:,:,:,:) = reshape(single(Z'),[DIM nScan]);
scans = arrayfun(@(i) sprintf('%s,%d',f,i),1:nScan,'UniformOutput',false)';

% parent model with SPM
pdir = fullfile(fxdir,'parent');
spec = [];
spec.dir = {pdir};
spec.timing = struct('units','secs','RT',TR,'fmri_t',16,'fmri_t0',8);
spec.sess = build_sess(scans,onsA,onsB,onsE,nuis);
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
mb = []; mb{1}.spm.stats.fmri_est.spmmat = {spmmat};
mb{1}.spm.stats.fmri_est.write_residuals = 0;
mb{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',mb);

% conditions of interest for LSS: A and B
SPM = load(spmmat).SPM;
cl = [];
for n = 1:2
    cl(n).sess = 1; cl(n).number = n; cl(n).name = char(SPM.Sess(1).U(n).name(1)); %#ok<AGROW>
    cl(n).file_name = ['[Run_1]_[Cond_' num2str(n) ']_[' regexprep(char(SPM.Sess(1).U(n).name(1)),' ','_') ']']; %#ok<AGROW>
end
fx = struct('maxmem',cs.maxmem,'spmmat',spmmat,'conditions',cl,'concat',false, ...
            'nscan',nScan,'session_mode','strict');
fprintf('fixture %s: %d in-brain voxels, %d interest trials, parent estimated: %s\n', ...
        cs.name,nnz(brain),nTrial,spmmat);
end

function s = build_sess(scans,onsA,onsB,onsE,R)
ons = {onsA,onsB,onsE}; nam = {'A','B','E'}; dur = {ones(numel(onsA),1),ones(numel(onsB),1),2*ones(numel(onsE),1)};
s.scans = scans;
for c = 1:3
    s.cond(c) = struct('name',nam{c},'onset',ons{c},'duration',dur{c}, ...
                       'tmod',0,'pmod',struct('name',{},'param',{},'poly',{}),'orth',1);
end
for q = 1:size(R,2), s.regress(q) = struct('name',sprintf('R%d',q),'val',R(:,q)); end
s.multi = {''}; s.multi_reg = {''}; s.hpf = 128;
end
