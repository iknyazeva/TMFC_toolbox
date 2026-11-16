function results = run_gppi_options_tests(out_dir,spm_dir,cases)

% Options agreement tests for the gPPI acceleration pooled-Cy / warm-start slice
% : run the CANDIDATE tmfc_gPPI / tmfc_gPPI_FIR three times on
% IDENTICAL fixture inputs (built once per case by the UNMODIFIED
% gppi_exact_fixture of tests/gppi/exact, called by name):
%   variant 1 - exact default (no option fields in tmfc.defaults)
%   variant 2 - gppi_method='pooled'
%   variant 3 - gppi_method='pooled' + gppi_reml_warm_start=true
% and compare 2-vs-1, 3-vs-1 and 3-vs-2 with the shared comparator
% (identical counts/names/shapes/NaN patterns are hard asserts; rel
% max/L2 are REPORT-ONLY: pooled/warm are approximations by design and the
% 1e-3 value is a target, not a gate -- evaluation/gppi/OPTIONS_IMPL.md).
% The exact-vs-main regression is run SEPARATELY by run_gppi_exact_tests.
%
% Pooled evidence per variant: the GLM_batches/pooled/<sub>_pooled_Cy.mat
% file must exist exactly in the pooled variants (2 and 3) and must not
% exist in the exact variant.
%
% FORMAT results = run_gppi_options_tests(out_dir,spm_dir,cases)
%   out_dir - disposable output folder (fresh per start)
%   spm_dir - SPM25 install
%   cases   - cell of exact-fixture case names (a2 cases; e.g.
%             {'ar1','fast','fir_ar1','fir_fast'})
%
% License: GPL-3.0-or-later.

if nargin < 3, error('Provide out_dir, spm_dir and cases.'); end
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));       % tests/gppi/options -> root
candidate = fullfile(repo,'functions');
exact_dir = fullfile(fileparts(here),'exact');
assert(exist(fullfile(candidate,'tmfc_gPPI.m'),'file') == 2, ...
    'candidate functions folder not found: %s',candidate);
assert(exist(fullfile(exact_dir,'gppi_exact_fixture.m'),'file') == 2, ...
    'exact fixture folder not found: %s',exact_dir);

oldpath = path; restore = onCleanup(@() path(oldpath));
addpath(spm_dir); addpath(exact_dir); addpath(here);
set(0,'DefaultFigureVisible','off');
spm('defaults','fmri'); spm_jobman('initcfg');
addpath(genpath(candidate),'-begin');
fixture_builder = which('gppi_exact_fixture');

opts = cell(3,1);
opts{1} = struct();                                          % exact default
opts{2} = struct('gppi_method','pooled');
opts{3} = struct('gppi_method','pooled','gppi_reml_warm_start',true);
vtag = {'exact','pooled','pooled_warm'};

if ~exist(out_dir,'dir'), mkdir(out_dir); end
results = struct([]);
for k = 1:numel(cases)
    name = cases{k};
    try
        [~,fx] = gppi_exact_fixture(name,fullfile(out_dir,'fixture',name));
        for v = 1:3
            proj = fullfile(out_dir,'runs',name,num2str(v));
            assert(~(exist(proj,'dir') == 7 && ~isempty(dir(fullfile(proj,'*')))), ...
                'Stale run tree %s from a previous attempt; use a fresh out_dir',proj);
            if ~exist(proj,'dir'), mkdir(proj); end
            copyfile(fullfile(fx.inputs_dir,'*'),proj);
            tmfc = gppi_options_tmfc(fx,proj,opts{v});
            assert(feval(fx.func,tmfc,1) == 1, ...
                'sub_check failed for %s variant %s',name,vtag{v});
        end

        r1 = fullfile(out_dir,'runs',name,'1');
        r2 = fullfile(out_dir,'runs',name,'2');
        r3 = fullfile(out_dir,'runs',name,'3');
        m21 = gppi_options_compare(r1,r2,fx);       % pooled vs exact (L2)
        m31 = gppi_options_compare(r1,r3,fx);       % warm vs exact (L3)
        m32 = gppi_options_compare(r2,r3,fx);       % warm vs pooled (L3)

        % pooled evidence files: present exactly in variants 2 and 3
        ev = false(1,3);
        for v = 1:3
            rv = fullfile(out_dir,'runs',name,num2str(v));
            ev(v) = exist(fullfile(rv,'ROI_sets',fx.set_name,fx.meth_dir, ...
                'GLM_batches','pooled',[fx.subject '_pooled_Cy.mat']),'file') == 2;
        end
        assert(~ev(1) && ev(2) && ev(3), ...
            'pooled Cy evidence file presence wrong for %s: %d%d%d', ...
            name,ev(1),ev(2),ev(3));

        e = struct('case_name',name,'method',fx.func,'noise',fx.cvi, ...
            'fixture_builder',fixture_builder,'error_msg','', ...
            'rel_max_pooled_vs_exact',m21.rel_max, ...
            'rel_l2_pooled_vs_exact',m21.rel_l2, ...
            'rel_max_warm_vs_exact',m31.rel_max, ...
            'rel_l2_warm_vs_exact',m31.rel_l2, ...
            'rel_max_warm_vs_pooled',m32.rel_max, ...
            'rel_l2_warm_vs_pooled',m32.rel_l2, ...
            'batch_equal_pooled_vs_exact',m21.batch_equal, ...
            'batch_equal_warm_vs_exact',m31.batch_equal, ...
            'batch_equal_warm_vs_pooled',m32.batch_equal, ...
            'n_matrix_finite',m21.n_matrix_finite, ...
            'pooled_cy_evidence',[ev(1) ev(2) ev(3)]);
        if isempty(results), results = e; else, results(end+1) = e; end %#ok<AGROW>
        save_progress(out_dir,results);
        fprintf(['%s (%s, %s): pooled-vs-exact rel_max=%g rel_l2=%g; ' ...
            'warm-vs-exact rel_max=%g rel_l2=%g; warm-vs-pooled rel_max=%g rel_l2=%g\n'], ...
            name,fx.func,fx.cvi,m21.rel_max,m21.rel_l2, ...
            m31.rel_max,m31.rel_l2,m32.rel_max,m32.rel_l2);
        if m21.rel_max > 1e-3 || m21.rel_l2 > 1e-3
            fprintf('  NOTE %s: pooled-vs-exact exceeds the 1e-3 REPORT-ONLY target\n',name);
        end
        if m31.rel_max > 1e-3 || m31.rel_l2 > 1e-3
            fprintf('  NOTE %s: warm-vs-exact exceeds the 1e-3 REPORT-ONLY target\n',name);
        end
    catch ME
        e = struct('case_name',name,'method','','noise','', ...
            'fixture_builder',fixture_builder, ...
            'error_msg',getReport(ME,'basic'));
        fields = {'rel_max_pooled_vs_exact','rel_l2_pooled_vs_exact', ...
            'rel_max_warm_vs_exact','rel_l2_warm_vs_exact', ...
            'rel_max_warm_vs_pooled','rel_l2_warm_vs_pooled', ...
            'batch_equal_pooled_vs_exact','batch_equal_warm_vs_exact', ...
            'batch_equal_warm_vs_pooled','n_matrix_finite','pooled_cy_evidence'};
        for f = fields
            e.(f{1}) = NaN; %#ok<AGROW>
        end
        if ~isfield(e,'pooled_cy_evidence') || numel(e.pooled_cy_evidence) ~= 3
            e.pooled_cy_evidence = [NaN NaN NaN];
        end
        if isempty(results), results = e; else, results(end+1) = e; end %#ok<AGROW>
        save_progress(out_dir,results);
        rethrow(ME);                     % stop; never skip a failure
    end
end
save_progress(out_dir,results);
end

function save_progress(out_dir,results)
save(fullfile(out_dir,'results.mat'),'results');
fid = fopen(fullfile(out_dir,'results.json'),'w');
fwrite(fid,jsonencode(results)); fclose(fid);
end
