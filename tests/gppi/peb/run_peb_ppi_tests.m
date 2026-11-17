function results = run_peb_ppi_tests(out_dir,ref_functions_dir,spm_dir,cases)

% Regression tests for the PEB deconvolution hoist in tmfc_PEB_PPI: compares the
% candidate (../../../functions) against a reference extraction of main
% (ref_functions_dir) on the real-SPM fixtures from peb_ppi_fixture.m.
%
% Gates per case:
%   - identical output file names and counts (== fixture expectation)
%   - PPI.name equality, PPI.ppi sizes equal, all values finite
%   - relative max |diff| AND relative L2 <= 1e-3 (no ULP requirement;
%     bitwise equality may hold in practice but is not asserted)
%   - spm_PEB call counts (measured via the MATLAB profiler FunctionTable,
%     documented fields FunctionName/NumCalls, the same method as the earlier
%     phase-timing harness) equal the fixture expectation: reference = nROI * (number
%     of selected conditions per run), candidate = nROI * (number of runs
%     with >= 1 selected condition) -- i.e. the expected reduction.
%     No path shadowing or spm_PEB stubs are used.
%
% FORMAT results = run_peb_ppi_tests(out_dir,ref_functions_dir,spm_dir,cases)
%   out_dir           - disposable output folder (fresh, outside sources)
%   ref_functions_dir - extraction of main:functions (read-only reference)
%   spm_dir           - SPM25 install (read-only)
%   cases             - optional subset of {'strict','repeated','centering_off'}
%
% STATUS: executed, PASS (see evaluation/gppi/REPORT.md, Agreement). Reproduction:
%
%   matlab -batch "cd <this_checkout>/tests/gppi/peb; results = ...
%       run_peb_ppi_tests('/tmp/tmfc10_peb_out', ...
%       '<main_toolbox>/functions', ...
%       '<spm25>'); disp(results)"
%
% License: GPL-3.0-or-later

if nargin < 3, error('Provide out_dir, ref_functions_dir and spm_dir.'); end
assert(exist(fullfile(ref_functions_dir,'tmfc_PEB_PPI.m'),'file') == 2);
assert(exist(fullfile(spm_dir,'spm_PEB.m'),'file') == 2);
here      = fileparts(mfilename('fullpath'));
candidate = fullfile(fileparts(fileparts(fileparts(here))),'functions');  % no '..' so which() prefix checks are exact
assert(exist(fullfile(candidate,'tmfc_PEB_PPI.m'),'file') == 2);
oldpath = path; restore = onCleanup(@() path(oldpath));
addpath(spm_dir); addpath(here); addpath(candidate);
spm('defaults','fmri'); spm_jobman('initcfg');
set(0,'DefaultFigureVisible','off');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if nargin < 4 || isempty(cases)
    cases = {'strict','repeated','centering_off'};
end

results = struct([]);
for k = 1:numel(cases)
    name = cases{k};
    [cs,fx] = peb_ppi_fixture(name,fullfile(out_dir,'fixture',name)); %#ok<ASGLU>

    % reference (main) and candidate (accelerated) runs, profiler off
    [ref_dir,w_ref]  = run_once(fx,ref_functions_dir,fullfile(out_dir,'runs',name,'ref'));
    [cand_dir,w_cand] = run_once(fx,candidate,        fullfile(out_dir,'runs',name,'cand'));

    % numeric gates
    m = compare_ppi_dirs(ref_dir,cand_dir);
    m.case_name = name;
    m.which_ref  = w_ref;
    m.which_cand = w_cand;
    assert(m.n_files == fx.expected_outputs,'Output count differs from fixture expectation');
    assert(m.finite,'Non-finite values in PPI outputs');
    assert(m.rel_max <= 1e-3 && m.rel_l2 <= 1e-3,'Relative max/L2 gate failed');

    % spm_PEB call-count gates (profiler FunctionTable, no shadowing)
    n_ref  = run_counted(fx,ref_functions_dir,fullfile(out_dir,'runs',name,'ref_prof'));
    n_cand = run_counted(fx,candidate,        fullfile(out_dir,'runs',name,'cand_prof'));
    assert(n_ref  == fx.expected_calls_old,'Reference spm_PEB call count unexpected');
    assert(n_cand == fx.expected_calls_new,'Candidate spm_PEB call count not reduced (PEB hoist)');
    m.calls_ref = n_ref; m.calls_cand = n_cand;

    if isempty(results), results = m; else, results(end+1) = m; end %#ok<AGROW>
    fprintf('%s: files=%d rel_max=%g rel_l2=%g finite=%d spm_PEB calls %d -> %d\n', ...
        name,m.n_files,m.rel_max,m.rel_l2,m.finite,n_ref,n_cand);
end
save(fullfile(out_dir,'results.mat'),'results');
end

%% ========================================================================

function [ppi_dir,used_file] = run_once(fx,func_dir,proj_dir)
% One tmfc_PEB_PPI run: writes identical VOI inputs, builds the CACHE the
% same way tmfc_PPI.m does, calls func_dir's tmfc_PEB_PPI, returns the
% output folder with the PPI_*.mat files and the resolved source file.
if exist(proj_dir,'dir'), rmdir(proj_dir,'s'); end
mkdir(proj_dir);
write_vois(fx,proj_dir);
cache = make_cache(fx,proj_dir);
addpath(func_dir,'-begin');
cleanup = onCleanup(@() rmpath(func_dir)); %#ok<NASGU>
clear tmfc_PEB_PPI
used_file = which('tmfc_PEB_PPI');
assert(~isempty(used_file) && strncmp(used_file,func_dir,numel(func_dir)), ...
    'tmfc_PEB_PPI did not resolve to %s (got %s)',func_dir,used_file);
fprintf('  [source] %s\n', used_file);
tmfc_PEB_PPI(cache,1);
ppi_dir = fullfile(proj_dir,'ROI_sets',fx.set_name,'PPIs','sub-01');
assert(exist(ppi_dir,'dir') == 7,'No PPI output directory produced');
end


function n = run_counted(fx,func_dir,proj_dir)
% One run under the profiler; returns the summed NumCalls of spm_PEB.
profile clear
profile on
cleanup = onCleanup(@() profile('off')); %#ok<NASGU>
run_once(fx,func_dir,proj_dir);
profile off
T = profile('info').FunctionTable;
hit = cellfun(@(nm) ~isempty(regexp(nm,'^spm_PEB$','once')),{T.FunctionName});
n   = sum([T(hit).NumCalls]);
end


function write_vois(fx,proj_dir)
% Per-run VOI files exactly where tmfc_PEB_PPI reads them. Old-style ROI 1
% has field Y only; 'ROI 2' also carries Yraw = Y + 1000, which must be
% ignored by tmfc_PEB_PPI (any accidental use fails the numeric gate).
voi_dir = fullfile(proj_dir,'ROI_sets',fx.set_name,'VOIs','sub-01');
if ~exist(voi_dir,'dir'), mkdir(voi_dir); end
for s = 1:size(fx.voi,1)
    for k = 1:numel(fx.roi_names)
        V = struct('Y',fx.voi{s,k});
        if ~fx.roi_old(k), V.Yraw = fx.voi{s,k} + 1000; end %#ok<STRNU>
        save(fullfile(voi_dir,sprintf('VOI_%s_%d.mat',fx.roi_names{k},s)),'-struct','V');
    end
end
end


function cache = make_cache(fx,proj_dir)
% CACHE built exactly as in tmfc_PPI.m (project_path, ROIset, conditions,
% session_mode, subjects); PPI_centering is fixed per fixture case.
cache = struct();
cache.project_path            = proj_dir;
cache.ROI_set_number          = 1;
cache.ROIset.set_name         = fx.set_name;
cache.ROIset.ROIs             = struct('name',{fx.roi_names{:}});
cache.ROIset.PPI_centering    = fx.centering;
cache.conditions.list         = fx.conditions;
cache.conditions.nCond        = numel(fx.conditions);
cache.session_mode            = fx.session_mode;
cache.subjects.sublist        = struct('path',fx.spmmat,'name','sub-01');
cache.subjects.nSub           = 1;
end


function m = compare_ppi_dirs(ref_dir,cand_dir)
a = dir(fullfile(ref_dir,'PPI_*.mat'));
b = dir(fullfile(cand_dir,'PPI_*.mat'));
assert(~isempty(a),'No reference outputs produced');
assert(isequal(sort({a.name}),sort({b.name})),'PPI file names differ');
rel_max = 0; rel_l2 = 0; finite_all = true;
for i = 1:numel(a)
    x = load(fullfile(a(i).folder,a(i).name));
    y = load(fullfile(b(i).folder,b(i).name));
    assert(strcmp(x.PPI.name,y.PPI.name),'PPI.name differs');
    assert(isequal(size(x.PPI.ppi),size(y.PPI.ppi)),'PPI.ppi size differs');
    finite_all = finite_all && all(isfinite(x.PPI.ppi(:))) && all(isfinite(y.PPI.ppi(:)));
    d = double(x.PPI.ppi(:)) - double(y.PPI.ppi(:));
    rel_max = max(rel_max, max(abs(d))/max(max(abs(double(x.PPI.ppi(:)))),realmin));
    rel_l2  = max(rel_l2,  norm(d)/max(norm(double(x.PPI.ppi(:))),realmin));
end
m = struct('n_files',numel(a),'rel_max',rel_max,'rel_l2',rel_l2,'finite',finite_all);
end
