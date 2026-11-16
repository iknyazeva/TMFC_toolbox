function results = run_lss_tests(out_dir,ref_functions_dir,spm_dir,cases)
% Regression against main. ref_functions_dir is a separate extraction of
% main:functions; out_dir is disposable. spm_dir identifies the SPM install.
% Tests exact float32 <= 1 ULP and max error / beta L2 <= 1e-3;
% rWLS exact uses the relative gate and reports float32 ULP separately;
% pooled uses correlations between ROI beta-series connectivity matrices.
% The caller must supply fresh paths outside the source tree.
if nargin < 3, error('Provide out_dir, ref_functions_dir and spm_dir.'); end
assert(exist(fullfile(ref_functions_dir,'tmfc_LSS.m'),'file') == 2);
here = fileparts(mfilename('fullpath'));
candidate = fullfile(here,'..','..','functions');
oldpath = path; restore = onCleanup(@() path(oldpath));
addpath(spm_dir); addpath(genpath(candidate)); addpath(here);
set(0,'DefaultFigureVisible','off');
spm('defaults','fmri'); spm_jobman('initcfg');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if nargin < 4 || isempty(cases)
    cases = {'none','ar1','fast','ar1_scaling_mask','single','collinear', ...
    'rank_def','concat','repeated','chunks','dense','var_std','var_dense', ...
    'rwls_single','rwls_outliers','rwls_concat'};
end
results = struct([]);
for k = 1:numel(cases)
    name = cases{k};
    if startsWith(name,'rwls_')
        [~,fx] = lss_rwls_fixture(name(6:end),fullfile(out_dir,'fixture',name));
        modes = {'exact'};
        if ~strcmp(name,'rwls_single'), modes = {'exact','pooled'}; end
    elseif startsWith(name,'var_')
        [~,fx] = lss_pooled_fixture(name,fullfile(out_dir,'fixture',name));
        modes = {'exact','pooled'};
    else
        [~,fx] = lss_exact_fixture(name,fullfile(out_dir,'fixture',name));
        modes = {'exact'};
    end
    checks = cell(1,numel(modes)+1);
    for i = 1:numel(checks)
        tmfc.project_path = fullfile(out_dir,'runs',name,num2str(i));
        if ~exist(tmfc.project_path,'dir'), mkdir(tmfc.project_path); end
        tmfc.subjects = struct('path',fx.spmmat,'name','sub-01');
        tmfc.defaults = struct('parallel',0,'maxmem',fx.maxmem,'resmem',true);
        tmfc.session_mode = fx.session_mode; tmfc.LSS.conditions = fx.conditions;
        if i == 1
            addpath(ref_functions_dir,'-begin');
        else
            addpath(candidate,'-begin');
            if ~strcmp(modes{i-1},'exact') || startsWith(name,'var_')
                tmfc.defaults.lss_method = modes{i-1};
            end
        end
        clear tmfc_LSS
        checks{i} = tmfc_LSS(tmfc);
        if i == 1, rmpath(ref_functions_dir); end
    end
    for i = 1:numel(modes)
        assert(isequaln(checks{1},checks{i+1}),'sub_check differs');
        ref = fullfile(out_dir,'runs',name,'1','LSS_regression','sub-01');
        got = fullfile(out_dir,'runs',name,num2str(i+1),'LSS_regression','sub-01');
        m = compare_outputs(ref,got);
        m.case_name = name; m.method = modes{i};
        if isempty(results), results = m; else, results(end+1) = m; end %#ok<AGROW>
        fprintf('%s %s: ulp=%g max/L2=%g ROI corr=%g\n',name,modes{i},m.ulp,m.relative,m.roi_corr);
        if strcmp(modes{i},'exact')
            assert(m.relative <= 1e-3,'Exact relative numerical gate failed');
            if ~startsWith(name,'rwls_')
                assert(m.ulp <= 1,'Exact float32 ULP gate failed');
            end
        else
            assert(isfinite(m.roi_corr) && m.roi_corr >= 0.999,'Pooled ROI correlation gate failed');
        end
    end
end
if any(strcmp(cases,'fast'))
    run_lss_pooled_options_tests(fullfile(out_dir,'pooled_options'),spm_dir);
end
save(fullfile(out_dir,'results.mat'),'results');
end

function m = compare_outputs(ref,got)
a = dir(fullfile(ref,'Betas','*.nii')); b = dir(fullfile(got,'Betas','*.nii'));
assert(~isempty(a) && isequal({a.name},{b.name}),'Beta file names differ');
ba = dir(fullfile(ref,'GLM_batches','*.mat')); bb = dir(fullfile(got,'GLM_batches','*.mat'));
assert(isequal({ba.name},{bb.name}) && numel(a)==numel(ba),'Batch file names differ');
for i = 1:numel(ba)
    x = load(fullfile(ba(i).folder,ba(i).name)); y = load(fullfile(bb(i).folder,bb(i).name));
    x.matlabbatch{1}.spm.stats.fmri_spec.dir = {''};
    y.matlabbatch{1}.spm.stats.fmri_spec.dir = {''};
    assert(isequaln(x,y),'Batch contents differ');
end
ulp = 0; relative = 0; ar = []; br = [];
for i = 1:numel(a)
    x = single(spm_read_vols(spm_vol(fullfile(a(i).folder,a(i).name))));
    y = single(spm_read_vols(spm_vol(fullfile(b(i).folder,b(i).name))));
    assert(isequal(isnan(x),isnan(y)),'Mask differs');
    valid = isfinite(x) & isfinite(y); assert(any(valid(:)),'Empty mask');
    assert(isequal(isfinite(x),isfinite(y)),'Finite masks differ');
    d = abs(double(x(valid))-double(y(valid)));
    relative = max(relative,max(d)/max(norm(double(x(valid))),realmin));
    ux = ordered(x(valid)); uy = ordered(y(valid));
    ulp = max(ulp,max(abs(double(ux)-double(uy))));
    % Eight fixed spatial ROIs; beta-series correlations define connectivity.
    [i1,i2,i3] = ind2sub(size(x),find(valid));
    labels = 1+(i1 > size(x,1)/2)+2*(i2 > size(x,2)/2)+4*(i3 > size(x,3)/2);
    xv = double(x(valid)); yv = double(y(valid));
    for j = 1:8
        assert(any(labels==j),'Empty ROI');
        ar(i,j) = mean(xv(labels==j)); br(i,j) = mean(yv(labels==j)); %#ok<AGROW>
    end
end
ca = corr(ar); cb = corr(br); edges = triu(true(8),1);
roi = corr(ca(edges),cb(edges));
m = struct('ulp',ulp,'relative',relative,'roi_corr',roi);
end

function u = ordered(x)
b = typecast(x(:),'uint32'); negative = bitget(b,32) ~= 0;
u = bitor(b,uint32(2147483648)); u(negative) = bitcmp(b(negative));
end
