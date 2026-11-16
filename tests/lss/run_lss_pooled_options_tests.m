function run_lss_pooled_options_tests(out_dir,spm_dir)
% Verify pooled defaults for FAST and AR(1), explicit overrides, and warnings.
if nargin < 2, error('Provide out_dir and spm_dir.'); end
here = fileparts(mfilename('fullpath')); addpath(here); addpath(spm_dir);
addpath(genpath(fullfile(here,'..','..','functions')));
spm('defaults','fmri'); spm_jobman('initcfg'); set(0,'DefaultFigureVisible','off');
for parent = {'fast','ar1'}
    name = parent{1}; root = fullfile(out_dir,name);
    [~,fx] = lss_exact_fixture(name,fullfile(root,'fixture'));
    modes = {'base','default_explicit','override','every4','exactwarn'};
    for k = 1:numel(modes)
        tmfc.project_path = fullfile(root,modes{k});
        if ~exist(tmfc.project_path,'dir'), mkdir(tmfc.project_path); end
        tmfc.subjects = struct('path',fx.spmmat,'name','sub-01');
        tmfc.defaults = struct('parallel',0,'maxmem',fx.maxmem,'resmem',true,'lss_method','pooled');
        tmfc.session_mode = fx.session_mode; tmfc.LSS.conditions = fx.conditions;
        switch modes{k}
            case 'default_explicit', tmfc.defaults.lss_reml_warm_start = strcmp(name,'fast');
            case 'override', tmfc.defaults.lss_reml_warm_start = ~strcmp(name,'fast');
            case 'every4', tmfc.defaults.lss_reml_warm_start = true; tmfc.defaults.lss_reml_every = 4;
            case 'exactwarn', tmfc.defaults.lss_method = 'exact'; tmfc.defaults.lss_reml_every = 4;
        end
        lastwarn('',''); tmfc_LSS(tmfc);
        if strcmp(modes{k},'exactwarn')
            [~,id] = lastwarn; assert(strcmp(id,'tmfc_LSS:PooledOptionsIgnored'),'Missing ignored-options warning');
        end
    end
    ref = fullfile(root,'base','LSS_regression','sub-01','Betas');
    a = dir(fullfile(ref,'*.nii'));
    for tag = {'default_explicit','override','every4'}
        b = dir(fullfile(root,tag{1},'LSS_regression','sub-01','Betas','*.nii'));
        assert(~isempty(a) && isequal({a.name},{b.name}),'File names differ');
        rel = 0; neq = 0;
        for k=1:numel(a)
            x=spm_read_vols(spm_vol(fullfile(a(k).folder,a(k).name)));
            y=spm_read_vols(spm_vol(fullfile(b(k).folder,b(k).name)));
            assert(isequal(isnan(x),isnan(y)),'Masks differ');
            q=isfinite(x)&isfinite(y); dx=double(x(q)); dy=double(y(q));
            rel=max(rel,max(abs(dx-dy))/max(norm(dx),realmin)); neq=neq+nnz(dx~=dy);
        end
        fprintf('OPTIONS %s %s max_rel %.12g unequal_voxels %d\n',name,tag{1},rel,neq);
        if strcmp(tag{1},'default_explicit'), assert(neq==0,'Pooled parent-specific default changed'); end
    end
end
end
