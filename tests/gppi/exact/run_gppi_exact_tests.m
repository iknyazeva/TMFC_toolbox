function results = run_gppi_exact_tests(out_dir,ref_functions_dir,spm_dir,cases)

% Exact-comparison tests for the gPPI batch/I/O changes (initcfg
% hoist + base_spec/VOI reuse) against main. Runs each fixture case twice
% over identical inputs: run 1 with the MAIN functions subtree on the path,
% run 2 with the CANDIDATE subtree; the baseline path (SPM + this folder)
% is restored between variants, so main never resolves candidate helpers
% (rwls helper isolation); changed function handles are cleared between
% variants. Parent fixtures are built on the candidate subtree (labeled in
% results as fixture_builder). Compares GLM batches (isequaln after root
% tokenization with preserved suffixes), ROI-to-ROI matrices (identical
% names/shape/NaN pattern, intentional NaN diagonals, finite off-diagonals)
% and seed-to-voxel beta images (identical names, geometry, NaN AND Inf
% masks). Numeric gate: relative max AND relative L2 <= 1e-3 (no ULP/
% bitwise gate). Progress is saved after EVERY case, before the next;
% failures are caught with the case identifier, recorded and rethrown
% (no skip-and-claim).
%
% FORMAT results = run_gppi_exact_tests(out_dir,ref_functions_dir,spm_dir,cases)
%   out_dir          - disposable output folder
%   ref_functions_dir- main checkout functions folder
%   spm_dir          - SPM install (SPM25)
%   cases            - optional cell array of fixture case names
%
% License: GPL-3.0-or-later

if nargin < 3, error('Provide out_dir, ref_functions_dir and spm_dir.'); end
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));        % tests/gppi/exact -> root
candidate = fullfile(repo,'functions');
assert(exist(fullfile(ref_functions_dir,'tmfc_gPPI.m'),'file') == 2, ...
    'ref_functions_dir must contain main tmfc_gPPI.m');
assert(exist(fullfile(candidate,'tmfc_gPPI.m'),'file') == 2, ...
    'candidate functions folder not found: %s',candidate);
oldpath = path; restore = onCleanup(@() path(oldpath));
addpath(spm_dir); addpath(here);
set(0,'DefaultFigureVisible','off');
spm('defaults','fmri'); spm_jobman('initcfg');
% baseline is captured AFTER initcfg so the matlabbatch/config directories
% added by SPM initialization stay on the path in every variant; variant
% isolation still holds (only SPM + this folder + the variant subtree)
baseline = path;
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if nargin < 4 || isempty(cases)
    cases = {'ar1','fir_ar1'};
end
results = struct([]);
for k = 1:numel(cases)
    name = cases{k};
    try
        % parent fixture built on the candidate subtree (input only)
        addpath(genpath(candidate),'-begin');
        fixture_builder = which('gppi_exact_fixture');
        [~,fx] = gppi_exact_fixture(name,fullfile(out_dir,'fixture',name));
        path(baseline);

        pending = false;   % this case IS executed in this invocation; the
        % not-yet-admitted scope lives in the reports, not in this field
        checks = cell(1,2); srcs = cell(1,2); hsrc = cell(1,2);
        for i = 1:2
            proj = fullfile(out_dir,'runs',name,num2str(i));
            % fresh-tree guard: leftovers from an earlier attempt would
            % inflate the comparator's dir() counts and fail confusingly;
            % use a fresh out_dir per start instead
            assert(~(exist(proj,'dir') == 7 && ~isempty(dir(fullfile(proj,'*')))), ...
                'Stale run tree %s from a previous attempt; use a fresh out_dir',proj);
            if ~exist(proj,'dir'), mkdir(proj); end
            copyfile(fullfile(fx.inputs_dir,'*'),proj);       % VOIs, PPIs, masks

            tmfc = struct();
            tmfc.project_path = proj;
            tmfc.subjects = struct('path',fx.spmmat,'name',fx.subject);
            tmfc.defaults = struct('parallel',0,'maxmem',fx.maxmem, ...
                'resmem',true,'analysis',fx.analysis);
            tmfc.session_mode = fx.session_mode;
            tmfc.ROI_set = struct('set_name',fx.set_name,'type','binary_images');
            for j = 1:numel(fx.roi_names)
                tmfc.ROI_set(1).ROIs(j).name = fx.roi_names{j};
                tmfc.ROI_set(1).ROIs(j).path_masked(1).subjects = ...
                    fullfile(proj,'masks',['ROI_' regexprep(fx.roi_names{j},' ','_') '_masked.nii']);
            end
            tmfc.ROI_set(1).gPPI.conditions = fx.conditions;
            if strcmp(fx.method,'fir')
                tmfc.ROI_set(1).gPPI_FIR.window = fx.window;
                tmfc.ROI_set(1).gPPI_FIR.bins   = fx.bins;
            end

            % Cheap existence checks before any MATLAB-heavy work
            assert(exist(fx.spmmat,'file') == 2,'parent SPM.mat missing');
            for j = 1:numel(fx.roi_names)
                assert(exist(tmfc.ROI_set(1).ROIs(j).path_masked(1).subjects,'file') == 2, ...
                    'copied mask missing for %s',fx.roi_names{j});
            end
            nv = numel(dir(fullfile(proj,'ROI_sets',fx.set_name,'VOIs',fx.subject,'VOI_*.mat')));
            npp = numel(dir(fullfile(proj,'ROI_sets',fx.set_name,'PPIs',fx.subject,'PPI_*.mat')));
            assert(nv == numel(fx.roi_names)*fx.nRun, ...
                'copied VOI count wrong: %d (expected %d)',nv,numel(fx.roi_names)*fx.nRun);
            assert(npp == numel(fx.roi_names)*numel(fx.conditions), ...
                'copied PPI count wrong: %d (expected %d)',npp,numel(fx.roi_names)*numel(fx.conditions));

            % variant path isolation: ONLY baseline + the variant subtree
            path(baseline);
            if i == 1
                addpath(genpath(ref_functions_dir),'-begin');
                sub_root = ref_functions_dir;
            else
                addpath(genpath(candidate),'-begin');
                sub_root = candidate;
            end
            clear tmfc_gPPI tmfc_gPPI_FIR tmfc_spm_rwls_spm ...
                tmfc_spm_rwls_reml tmfc_spm_rwls_est_non_sphericity
            sub_root_sep = [sub_root filesep];   % boundary: a sibling dir
            % extending sub_root (e.g. <root>_backup) must NOT pass
            srcs{i} = which(fx.func);
            assert(~isempty(srcs{i}) && strncmp(srcs{i},sub_root_sep,numel(sub_root_sep)), ...
                'run %d resolved %s from %s (expected under %s)', ...
                i,fx.func,srcs{i},sub_root);
            h = helper_paths();
            if strcmp(fx.cvi,'wls')
                assert(~isempty(h.rwls_spm) && strncmp(h.rwls_spm,sub_root_sep,numel(sub_root_sep)), ...
                    'RWLS helper not under the variant root for %s run %d: %s', ...
                    name,i,h.rwls_spm);
            end
            hsrc{i} = sprintf('rwls_spm=%s;rwls_reml=%s;rwls_est_ns=%s;est_ns=%s', ...
                h.rwls_spm,h.rwls_reml,h.rwls_est_ns,h.est_ns);
            checks{i} = feval(fx.func,tmfc,1);
            assert(checks{i} == 1,'sub_check failed for %s run %d',name,i);
        end
        path(baseline);

        ref = fullfile(out_dir,'runs',name,'1');
        got = fullfile(out_dir,'runs',name,'2');
        m = gppi_exact_compare(ref,got,fx);   % SAME comparator as verify_gppi_exact_saved
        m.case_name = name; m.method = fx.func; m.pending = pending; %#ok<STRUCT>
        m.noise = fx.cvi; m.analysis = fx.analysis; %#ok<STRUCT>
        m.ref_source = srcs{1}; m.candidate_source = srcs{2}; %#ok<STRUCT>
        m.helper_ref_source = hsrc{1}; m.helper_candidate_source = hsrc{2}; %#ok<STRUCT>
        m.fixture_builder = fixture_builder; m.error_msg = ''; %#ok<STRUCT>
        if isempty(results), results = m; else, results(end+1) = m; end %#ok<AGROW>
        save_progress(out_dir,results);                % progress BEFORE next case
        if pending, pend_s = ' PENDING'; else, pend_s = ''; end
        fprintf('%s: ref=%s cand=%s\n  helpers ref: %s\n  helpers cand: %s\n', ...
            name,srcs{1},srcs{2},hsrc{1},hsrc{2});
        fprintf('%s (%s%s): rel_max=%g rel_l2=%g batch_eq=%d\n',name,fx.func, ...
            pend_s,m.rel_max,m.rel_l2,m.batch_equal);
        assert(m.batch_equal,'GLM batch contents differ');
        assert(m.rel_max <= 1e-3 && m.rel_l2 <= 1e-3,'Numeric gate failed');
    catch ME
        if ~isempty(results) && strcmp(results(end).case_name,name)
            results(end).error_msg = getReport(ME,'basic');
        else
            e = result_template(name); e.error_msg = getReport(ME,'basic');
            if isempty(results), results = e; else, results(end+1) = e; end %#ok<AGROW>
        end
        save_progress(out_dir,results);
        rethrow(ME);                                  % stop; never skip a failure
    end
end
save_progress(out_dir,results);
end

function save_progress(out_dir,results)
save(fullfile(out_dir,'results.mat'),'results');
fid = fopen(fullfile(out_dir,'results.json'),'w');
fwrite(fid,jsonencode(results)); fclose(fid);
end

function m = result_template(name)
m = struct('rel_max',NaN,'rel_l2',NaN,'batch_equal',false, ...
    'n_matrix_finite',0,'n_beta_voxels',0,'case_name',name,'method','', ...
    'pending',true,'noise','','analysis',NaN,'ref_source','', ...
    'candidate_source','','helper_ref_source','', ...
    'helper_candidate_source','','fixture_builder','','error_msg','');
end

function h = helper_paths()
h.rwls_spm    = which('tmfc_spm_rwls_spm');
h.rwls_reml   = which('tmfc_spm_rwls_reml');
h.rwls_est_ns = which('tmfc_spm_rwls_est_non_sphericity');
h.est_ns      = which('spm_est_non_sphericity');
end
