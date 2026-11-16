# Reproducing the evaluation

## Inputs

For the paper simulations, use the `TMFC_simulations` time-series files,
task-design files and image-generation helpers, together with the relevant
`TMFC_SCRIPTS` pipeline. The table uses subject 1 and five designs:

| Design | ISI, s | Trials per condition | Seed id |
|---|---:|---:|---:|
| 02 | 6 | 100 | 2 |
| 02v, variable HRF | 6 | 100 | 202 |
| 12 | 2 | 100 | 12 |
| 15 | 12 | 100 | 15 |
| 17 | 6 | 40 | 17 |

TR is 2 s, event duration 1 s. For each design use
`rng(1000*design_id + subject)`. Drop the first three dummy scans, z-score
the subject's oscillation and coactivation series separately, add 100 to the
oscillation, divide coactivation by 1, and add `randn / 0.4`. Generate the
100-ROI images and shift task onsets by −0.2 s. Estimate a classical first
level with canonical HRF, no derivatives, AR(1), HPF 128 s, UFp 1, global
`'None'`, no explicit mask, TR 2 s, 16 microtime bins and reference bin 8.
The 02v input uses the variable-HRF simulation file with the 02 task design.

For ds000008, use sub-01 Stop Signal Task runs 1–3: 182 scans per run, TR 2 s,
4 mm FWHM SPM-smoothed fMRIPrep-derived BOLD in MNI space at 2 mm resolution,
and the corresponding first-level design and confounds. The model has 284
`go_correct`, 52 `stop_success`, and 44 `stop_failed` LSS trials, a nuisance
condition and 31–34 nuisance regressors per run. The parent uses HPF 128 s,
canonical HRF and UFp 0.001; its estimated mask contains 219,177 voxels.
Use parent variants for AR(1), independent errors, FAST and rWLS. The BSC
comparison uses nine available masks in the 10ROIs set and first eigenvariates.
The local smooth4 inputs have not been proven identical to another site's data.

The recorded measurements used a MacBook Pro M5 Pro (18 cores, 24 GB),
MATLAB R2025a and SPM25, with serial LSS and default BLAS threads. Synthetic
subject-1 elapsed time and elementwise maximum absolute z difference were
not recorded; they are marked not measured. Data and binary results are not
included here.

## Run the comparisons

Prepare separate writable project copies for `main`, `accelerated_lss` exact
and pooled, and `main_global_fixed` when comparing after-FIR results. Use the
same inputs and parent models in each copy. Set `tmfc.project_path` and the
subject paths for each copy; LSS writes beta images to the project directory.
Put SPM and the appropriate toolbox branch on the MATLAB path. Load your
project's `tmfc` configuration using your normal workflow, then run:

```matlab
% In a session using main:
tmfc.defaults.parallel = 0;
tmfc_LSS(tmfc);

% In a separate session using accelerated_lss and another project copy:
tmfc.defaults.parallel = 0;
tmfc.defaults.lss_method = 'exact';
tmfc_LSS(tmfc);

% In another copy with the same parents:
tmfc.defaults.lss_method = 'pooled';
tmfc_LSS(tmfc);
```

For FAST, pooled enables warm start by default; set
`tmfc.defaults.lss_reml_warm_start = false` to compare without it. Run
`tmfc_LSS_after_FIR(tmfc)` on corresponding project copies for after-FIR
comparisons. For rWLS, `exact` and `pooled` are supported; warm start and the
interval are ignored for rWLS with a warning. Compare corresponding
beta files and BSC results, retaining their masks and trial ordering.
Float32 outputs cannot establish a double-precision gate.

Run the fixture suite with a fresh output directory outside the source tree,
a separate extraction of the `main` functions, and the SPM directory:

```matlab
addpath('<toolbox>/tests/lss');
results = run_lss_tests('<disposable-output>', '<main-functions>', '<spm>');
```

 Timing, peak RSS and detailed agreement tables are in
[sub-01](sub-01/); the experimental FAST measurements are in
[fast_speedup.md](fast_speedup.md). The component benchmarks and the end-to-end
checks have different measurement protocols, as indicated in the
tables. Reproduction commands have not been rerun for this document.
