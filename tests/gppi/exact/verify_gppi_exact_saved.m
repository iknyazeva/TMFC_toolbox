function results = verify_gppi_exact_saved(out_root,cases)

% Verify ALREADY SAVED run_gppi_exact_tests outputs with the SAME
% comparator (gppi_exact_compare) used by the live runner. Executes NO
% GLM/estimation, regenerates nothing, writes nothing into the saved run
% trees: for each case it compares out_root/runs/<case>/1 (ref = main)
% vs out_root/runs/<case>/2 (candidate) and applies the live gates.
%
% Intended continuation after an interrupted run (e.g. the exact-test output
% final1: ar1 outputs are complete and immutable, the runner crashed in
% its REPORTING line before the gates): this validates them without any
% rerun of tmfc_gPPI/tmfc_gPPI_FIR.
%
% FORMAT results = verify_gppi_exact_saved(out_root,cases)
%   out_root - folder containing runs/<case>/{1,2}
%   cases    - optional cell array (default {'ar1','fir_ar1'}; missing
%              run trees for a case are an error)
%
% Saves out_root/results_verify.{mat,json}. Gates: batch equality after
% project-root tokenization with preserved suffixes; identical expected
% counts, names, shapes, NaN patterns, intentional NaN diagonals;
% rel max AND rel L2 <= 1e-3 (no ULP/bitwise gate).
%
% STATUS: executed, PASS (evaluation/gppi/GPPI_EXACT_TESTS.md).
%
% License: GPL-3.0-or-later

if nargin < 2 || isempty(cases), cases = {'ar1','fir_ar1'}; end
results = struct([]);
for k = 1:numel(cases)
    name = cases{k};
    ref = fullfile(out_root,'runs',name,'1');
    got = fullfile(out_root,'runs',name,'2');
    assert(exist(ref,'dir') == 7 && exist(got,'dir') == 7, ...
        'saved run trees missing for %s',name);
    fx = gppi_exact_case_meta(name);
    m = gppi_exact_compare(ref,got,fx);
    m.case_name = name; m.method = fx.func; m.noise = fx.cvi; %#ok<STRUCT>
    m.analysis = fx.analysis; m.pending = false; %#ok<STRUCT>
    if isempty(results), results = m; else, results(end+1) = m; end %#ok<AGROW>
    fprintf('%s (%s): rel_max=%g rel_l2=%g batch_eq=%d\n', ...
        name,fx.func,m.rel_max,m.rel_l2,m.batch_equal);
    assert(m.batch_equal,'GLM batch contents differ');
    assert(m.rel_max <= 1e-3 && m.rel_l2 <= 1e-3,'Numeric gate failed');
end
save(fullfile(out_root,'results_verify.mat'),'results');
fid = fopen(fullfile(out_root,'results_verify.json'),'w');
fwrite(fid,jsonencode(results)); fclose(fid);
end
