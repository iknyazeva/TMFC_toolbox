function m = gppi_exact_compare(ref,got,fx)

% Shared comparator for gPPI acceleration gPPI exact tests: the SAME
% comparator is used by the live runner (run_gppi_exact_tests) and by the
% verify-only path (verify_gppi_exact_saved), so already saved outputs can
% be validated later with identical logic. No GLM/estimation is executed.
%
% FORMAT m = gppi_exact_compare(ref,got,fx)
%   ref, got - project roots of the two runs (ref = main, got = candidate)
%   fx       - case metadata bundle (gppi_exact_fixture at run time, or
%              gppi_exact_case_meta for saved outputs): set_name, subject,
%              roi_names, meth_dir, analysis, expected_batches,
%              expected_matrices, expected_betas_per_roi
%   m        - rel_max, rel_l2, batch_equal, n_matrix_finite, n_beta_voxels
%
% Gates applied by the CALLER: batch_equal; rel max AND rel L2 <= 1e-3.
% Here: expected counts, identical names; matrices identical shape/NaN
% pattern with NaN ONLY on the diagonal (by design); beta images identical
% NaN masks; GLM batches isequaln after tokenizing ONLY the known
% per-run project-root prefix with the SAME token, PRESERVING the
% relative suffix (mask/filename differences still fail).
%
% License: GPL-3.0-or-later

base = @(proj) fullfile(proj,'ROI_sets',fx.set_name,fx.meth_dir);
rel_max = 0; rel_l2 = 0; batch_equal = true; n_matrix_finite = 0; n_beta_voxels = 0;

% Roots must be disjoint (neither a path-prefix of the other): otherwise
% the longer root would be tokenized against the shorter one, producing a
% WRONG relative suffix. Per-run project roots are siblings by design
% (<out>/runs/<case>/{1,2}); this guard makes misuse loud, not silent.
assert(isequal(ref,got) || ~strncmp(ref,got,min(numel(ref),numel(got))), ...
    'Comparator roots must be disjoint siblings: %s vs %s',ref,got);

% GLM batches: expected count, identical names and contents after
% normalizing known output-root paths (fmri_spec.dir and any field that
% carries the per-run project root); all other fields must match
ba = dir(fullfile(base(ref),'GLM_batches','*',[fx.subject '_*_GLM.mat']));
bb = dir(fullfile(base(got),'GLM_batches','*',[fx.subject '_*_GLM.mat']));
assert(numel(ba) == fx.expected_batches && numel(bb) == fx.expected_batches, ...
    'Batch count differs: %d vs %d (expected %d)',numel(ba),numel(bb),fx.expected_batches);
assert(isequal({ba.name},{bb.name}),'Batch file names differ');
for i = 1:numel(ba)
    x = load(fullfile(ba(i).folder,ba(i).name));
    y = load(fullfile(bb(i).folder,bb(i).name));
    xn = normalize_roots(x.matlabbatch,ref,got);
    yn = normalize_roots(y.matlabbatch,ref,got);
    batch_equal = batch_equal && isequaln(xn,yn);
end

% ROI-to-ROI matrices (analysis 1 or 2): expected count, identical
% names/shape/NaN pattern; NaN only on the diagonal (by design)
if fx.analysis == 1 || fx.analysis == 2
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
end

% Seed-to-voxel beta images (analysis 1 or 3): expected count per ROI,
% identical names/masks
if fx.analysis == 1 || fx.analysis == 3
    for j = 1:numel(fx.roi_names)
        roi = fx.roi_names{j};
        na = dir(fullfile(base(ref),'Seed_to_voxel',roi,'*.nii'));
        nb = dir(fullfile(base(got),'Seed_to_voxel',roi,'*.nii'));
        assert(numel(na) == fx.expected_betas_per_roi && ...
            numel(nb) == fx.expected_betas_per_roi, ...
            'Beta count differs for %s: %d vs %d (expected %d)', ...
            roi,numel(na),numel(nb),fx.expected_betas_per_roi);
        assert(isequal({na.name},{nb.name}),'Beta file names differ');
        for i = 1:numel(na)
            va = spm_vol(fullfile(na(i).folder,na(i).name));
            vb = spm_vol(fullfile(nb(i).folder,nb(i).name));
            assert(isequal([va.dim],[vb.dim]) && isequal(va(1).mat,vb(1).mat), ...
                'Beta geometry differs');
            x = single(spm_read_vols(va));
            y = single(spm_read_vols(vb));
            assert(isequal(size(x),size(y)),'Beta shape differs');
            assert(isequal(isnan(x),isnan(y)),'Beta NaN mask differs');
            assert(isequal(isinf(x),isinf(y)),'Beta Inf mask differs');
            v = ~(isnan(x)|isinf(x));
            assert(any(v(:)),'Empty finite beta mask');
            n_beta_voxels = n_beta_voxels + nnz(v);
            d = abs(double(x(v))-double(y(v)));
            rel_max = max(rel_max,max(d)/max(max(abs(double(x(v)))),realmin));
            rel_l2  = max(rel_l2,norm(d)/max(norm(double(x(v))),realmin));
        end
    end
end
m = struct('rel_max',rel_max,'rel_l2',rel_l2,'batch_equal',batch_equal, ...
    'n_matrix_finite',n_matrix_finite,'n_beta_voxels',n_beta_voxels);
end

function out = normalize_roots(v,ref,got)
% Recursively tokenize ONLY the known per-run project-root prefix of char
% fields (output dirs, copied PPI/VOI/mask paths) with the SAME token,
% PRESERVING the relative suffix so mask/filename differences still fail;
% everything else is compared verbatim
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
