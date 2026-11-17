# Reproducing the gPPI evaluation

## Inputs

For the synthetic part, use the `TMFC_simulations` time-series files,
task-design files and image-generation helpers. The table uses subject 1
and four designs:

| Design | ISI, s | Trials per condition | Seed id |
|---|---|---:|---:|
| 01 block | 20 s blocks, 10 per condition | 10 | 1 |
| 01v block, variable HRF | 20 s blocks, 10 per condition | 10 | 201 |
| 12 event | 2 | 100 | 12 |
| 13 event | 4 | 100 | 13 |

TR is 2 s (event duration 1 s for 12/13). For each design use
`rng(1000*seed_id + subject)`. Drop the first three dummy scans, z-score
the subject's oscillation and coactivation series separately, add 100 to
the oscillation, divide coactivation by 1, and add `randn / 0.4`. Generate
the 100-ROI single-voxel images and shift task onsets by -0.2 s. For 01v
use the variable-HRF time-series file (per-subject coactivation tensor)
with the 01 task design. Estimate a classical first level with canonical
HRF, no derivatives, AR(1), HPF 128 s, UFp 1, global 'None', no explicit
mask, 16 microtime bins and reference bin 8. The ground truth for the A-B
contrast is `ground_truth_symm_matrix.mat` in the simulations repository.

For the real-data part, use ds000008 sub-01 Stop Signal Task runs 1-3
(182 scans per run, TR 2 s) with the 10ROIs masks and the precomputed
first-level parent as in the LSS evaluation; parent variants for AR(1)
and FAST are used here.

Measurements were recorded with MATLAB R2025a and SPM25, serial runs
without parallel workers, on one machine; each timing cell is a single
run. Data and binary results are not included.

## Run the comparisons

Put SPM and the toolbox on the MATLAB path. The tests build their own
small fixtures; the runners take explicit paths for the reference
(main) and candidate (this) `functions/` folders:

```matlab
% exact equivalence against main, real-SPM fixtures (14 executable cases)
run_gppi_exact_tests('<out>', '<main_toolbox>/functions', '<spm25>', ...
    {'ar1','fir_ar1','fast','fir_fast','rwls','fir_rwls', ...
     'ar1_a1','fir_ar1_a1','ar1_s2v','fir_ar1_s2v', ...
     'none_s2v','none_a1','fir_none_s2v','fir_none_a1'})

% compact option-guard preflight (exact-only build)
ok = test_gppi_option_guards('<this_checkout>/functions');

% PEB deconvolution hoist
run_peb_ppi_tests(...)        % see tests/gppi/peb/README.md
```

For the timing and pooled-agreement measurements, prepare two writable
project copies (main and this branch) with identical inputs and parent
models, and time the VOI, PPI, gPPI and gPPI-FIR branches separately (the
recorded numbers report branches, never sums). Historical pooled timings
were taken on pre-removal revisions with the option present.

## Recorded results

- `REPORT.md` — summary, parameters, guards, pooled math, timing and
  agreement tables (the colleague-facing overview).
- `sub-01/` — real-data timings and agreement per noise model.
- `simulations/sim_sub01.md` / `.csv` — synthetic single-subject table
  (AUC, agreement, per-branch timings).
