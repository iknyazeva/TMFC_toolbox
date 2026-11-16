function m = gppi_options_compare(ref,got,fx)

% Comparator for the gPPI acceleration options tests : compares two CANDIDATE
% run trees (exact vs pooled vs pooled+warm) of one analysis-2 fixture
% case. Adapted from tests/gppi/exact/gppi_exact_compare.m (NOT modified):
% the same GLM-batch equality after per-run project-root tokenization and
% the same ROI-to-ROI matrix checks (identical counts/names/shapes/NaN
% patterns, NaN only on the diagonal, finite off-diagonals). The numeric
% rel max/L2 are REPORT-ONLY here (approximation by design): the caller
% records them and only flags the 1e-3 target, it does not gate on them.
%
% FORMAT m = gppi_options_compare(ref,got,fx)
%   ref, got - project roots of the two candidate runs (sibling trees)
%   fx       - fixture bundle (gppi_exact_fixture): set_name, subject,
%              roi_names, meth_dir, expected_batches, expected_matrices
%   m        - rel_max, rel_l2, batch_equal, n_matrix_finite
%
% License: GPL-3.0-or-later.

assert(fx.analysis == 2,'gppi_options_compare: analysis-2 cases only');
base = @(proj) fullfile(proj,'ROI_sets',fx.set_name,fx.meth_dir);
rel_max = 0; rel_l2 = 0; batch_equal = true; n_matrix_finite = 0;

% Roots must be disjoint siblings (same guard as gppi_exact_compare)
assert(isequal(ref,got) || ~strncmp(ref,got,min(numel(ref),numel(got))), ...
    'Comparator roots must be disjoint siblings: %s vs %s',ref,got);

% GLM batches: identical count/names/contents after root tokenization
ba = dir(fullfile(base(ref),'GLM_batches','*',[fx.subject '_*_GLM.mat']));
bb = dir(fullfile(base(got),'GLM_batches','*',[fx.subject '_*_GLM.mat']));
assert(numel(ba) == fx.expected_batches && numel(bb) == fx.expected_batches, ...
    'Batch count differs: %d vs %d (expected %d)', ...
    numel(ba),numel(bb),fx.expected_batches);
assert(isequal({ba.name},{bb.name}),'Batch file names differ');
for i = 1:numel(ba)
    x = load(fullfile(ba(i).folder,ba(i).name));
    y = load(fullfile(bb(i).folder,bb(i).name));
    xn = normalize_roots(x.matlabbatch,ref,got);
    yn = normalize_roots(y.matlabbatch,ref,got);
    batch_equal = batch_equal && isequaln(xn,yn);
end

% ROI-to-ROI matrices (asym + symm), NaN diagonal by design
fields = {'asymmetrical','ppi_matrix';'symmetrical','symm_ppi_matrix'};
for sd = 1:size(fields,1)
    ma = dir(fullfile(base(ref),'ROI_to_ROI',fields{sd,1},[fx.subject '_*.mat']));
    mb = dir(fullfile(base(got),'ROI_to_ROI',fields{sd,1},[fx.subject '_*.mat']));
    assert(numel(ma) == fx.expected_matrices && numel(mb) == fx.expected_matrices, ...
        'Matrix count differs in %s: %d vs %d (expected %d)', ...
        fields{sd,1},numel(ma),numel(mb),fx.expected_matrices);
    assert(isequal({ma.name},{mb.name}),'Matrix file names differ');
    for i = 1:numel(ma)
        xa = load(fullfile(ma(i).folder,ma(i).name));
        xb = load(fullfile(mb(i).folder,mb(i).name));
        A = xa.(fields{sd,2}); B = xb.(fields{sd,2});
        assert(isequal(size(A),size(B)),'Matrix shape differs');
        assert(isequal(isnan(A),isnan(B)),'Matrix NaN pattern differs');
        assert(all(isnan(diag(A))) && all(isfinite(A(~logical(eye(size(A))))), 'all'), ...
            'NaN expected only on the matrix diagonal');
        v = ~isnan(A);
        n_matrix_finite = n_matrix_finite + nnz(v);
        d = abs(A(v)-B(v));
        rel_max = max(rel_max,max(d)/max(max(abs(A(v))),realmin));
        rel_l2  = max(rel_l2,norm(d)/max(norm(A(v)),realmin));
    end
end
m = struct('rel_max',rel_max,'rel_l2',rel_l2,'batch_equal',batch_equal, ...
    'n_matrix_finite',n_matrix_finite);
end

function out = normalize_roots(v,ref,got)
% Local copy of the normalize_roots subfunction of
% tests/gppi/exact/gppi_exact_compare.m (not modified): tokenize ONLY the
% known per-run project-root prefix of char fields with the SAME token,
% preserving the relative suffix.
if isstruct(v)
    out = v;
    for f = fieldnames(v)'
        for k = 1:numel(v)
            out(k).(f{1}) = normalize_roots(v(k).(f{1}),ref,got);
        end
    end
elseif iscell(v)
    out = cell(size(v));
    for k = 1:numel(v)
        out{k} = normalize_roots(v{k},ref,got); %#ok<AGROW>
    end
elseif ischar(v)
    out = v;
    if strncmp(v,ref,numel(ref))
        out = ['<ROOT>' v(numel(ref)+1:end)];
    elseif strncmp(v,got,numel(got))
        out = ['<ROOT>' v(numel(got)+1:end)];
    end
else
    out = v;
end
end
