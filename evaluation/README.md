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
[sub-01](lss/sub-01/); the experimental FAST measurements are in
[fast_speedup.md](lss/fast_speedup.md). The component benchmarks and the end-to-end
checks have different measurement protocols, as indicated in the
tables. Reproduction commands have not been rerun for this document.

## Evidence and materials

- LSS: [sub-01 evidence](lss/sub-01/), [synthetic simulations](lss/simulations/),
  [FAST speedup notes](lss/fast_speedup.md).
- gPPI: [ds000008 sub-01 timing and agreement](gppi/sub-01/),
  [synthetic single subject](gppi/simulations/).

## gPPI acceleration

### What you need (gPPI branch)

- MATLAB (R2025a was used) and SPM25 on the MATLAB path.
- This toolbox branch on the path (the `functions/` folder of this
  checkout).
- The `main` reference branch (`functions/` folder) for equivalence
  comparisons.

All tests below write only to an output directory you choose. They run in
a single MATLAB process without parallel workers. Output directories must
not exist beforehand.

### Fixture-based equivalence tests (no external data needed)

The fixtures are self-contained: `gppi_exact_fixture.m` builds a small
simulated dataset with a real SPM25 first-level model (AR(1)/FAST/rWLS/
independent-error parents) and precomputed PPI regressors. The runner
compares the candidate against the `main` reference on identical inputs.

```matlab
% add this checkout and the reference to the path
addpath('<this_checkout>/functions');
addpath('<this_checkout>/tests/gppi/exact');

% 22-case validation (all noise models x analyses x methods, except
%   analysis 2 with independent-error parents, which is excluded here;
%   its identical legacy failure is separate evidence, not one of
%   these 22 cases)
r = run_gppi_exact_tests('<fresh_out_dir>', ...
    '<main_checkout>/functions', ...          % reference
    '<spm25>', ...                             % SPM25 root
    {'ar1','fir_ar1','fast','fir_fast','rwls','fir_rwls', ...
     'ar1_a1','fir_ar1_a1','ar1_s2v','fir_ar1_s2v', ...
     'none_s2v','none_a1','fir_none_s2v','fir_none_a1', ...
     'fast_s2v','fast_a1','rwls_s2v','rwls_a1', ...
     'fir_fast_s2v','fir_fast_a1','fir_rwls_s2v','fir_rwls_a1'});

% option-guard preflight (exact-only delivery; ~seconds)
test_gppi_option_guards('<this_checkout>/functions');

% PEB deconvolution hoist (3 fixtures; ~1 min)
addpath('<this_checkout>/tests/gppi/peb');
run_peb_ppi_tests('<fresh_peb_out>', ...
    '<main_checkout>/functions', '<spm25>');
```

The full suite -- the 22 cases, the option-guard preflight and the three
PEB fixtures -- ran in one serial process in about 14 minutes; the
22-case block alone was not timed separately.

Results are read from the saved `<out_dir>/results.json` (per-case
rel max / rel L2, batch equality, structural checks, per-variant
`which` source attribution) and `<peb_out>/results.mat`. The pass gate
is relative max AND L2 <= 1e-3 on every case, with identical file
counts, names, shapes and NaN patterns.

The `gppi_exact_fixture.m` case list and seeds:
`ar1` (default, seed 11), s2v cases seed 12, a1 cases seed 13;
`none` uses independent-error parents, `fast` the FAST noise model,
`rwls` real robust estimation via `tmfc_spm_rwls_spm` (one covariance
component per scan; the parent is re-estimated in a clean subdirectory
to avoid the interactive overwrite guard). Old-style VOI files
(`Y` only, no `Yraw`) are covered by one ROI in every fixture.

### Real-data timing (ds000008 sub-01)

The recorded AR(1) timings ran the `main` branch first and the
accelerated branch second, sequentially in ONE `matlab -batch` process
(a launch driver), each branch on its own writable project copy, with
identical inputs, masks and parent models in the two copies (parent SPM
at AR(1) or FAST, 9 ROIs of
the `10ROIs` set, 9 conditions = 3 condition types × 3 runs (strict
session numbering, per the recorded manifest), ROI-to-ROI analysis 2,
strict session mode). Because the two runs shared one process, the
second run may have benefited from cache/JIT warm-up; the record does
not quantify this. The process layout of the recorded FAST timings is
not documented in the available records, so the one-process statement
is made for AR(1) only. Each branch was timed with the MATLAB profiler per phase (VOI,
PPI, gPPI, gPPI-FIR separately; the reports never sum phases). The
harness that produced the recorded numbers -- the two-subtree isolation,
phase timing and output comparison -- is NOT part of this repository (a
per-launch wrapper with a watchdog, snapshots and ledgers); the
commands below give an equivalent measurement procedure, not a
reproduction of the historical times. To measure:

1. Prepare two writable copies of your project with identical inputs,
   masks and parent models, one per branch, and configure each copy.
2. Put SPM25 and the appropriate `functions/` subtree on the path and
   load that copy's `tmfc` configuration. To match the recorded AR(1)
   protocol, run `main` first and the accelerated branch second,
   sequentially in one MATLAB process. For cleaner isolation, prefer a
   separate fresh MATLAB session per branch -- a different protocol,
   not identical to the recorded one.
3. Time each phase separately with `profile on; ...; profile off` and
   read `profile('info')`, or use `tic`/`toc` around each call:

```matlab
% Per branch: SPM25 and that branch's functions/ on the path, tmfc
% loaded for that branch's project copy. Recorded AR(1) protocol: both
% branches sequentially in ONE process (main first); a fresh session
% per branch is the cleaner alternative, not an identical protocol:
tmfc.defaults.parallel = 0;                % serial, as recorded
t = tic; tmfc_VOI(tmfc);      voi_s  = toc(t);
t = tic; tmfc_PPI(tmfc);      ppi_s  = toc(t);
t = tic; tmfc_gPPI(tmfc);     gppi_s = toc(t);
t = tic; tmfc_gPPI_FIR(tmfc); fir_s  = toc(t);
```

4. Compare the outputs of the two runs per saved matrix: check identical
   file names, shapes and NaN masks, then compute rel max / rel L2 over
   the finite off-diagonal entries, candidate `a` vs reference `b`, as
   `max(abs(a-b))/max(max(abs(b)),realmin)` and
   `norm(a-b)/max(norm(b),realmin)` (`realmin` guards the zero
   denominator). The shipped fixture comparator
   `tests/gppi/exact/gppi_exact_compare.m` expects the fixture output
   layout (expected batch/beta structure) and is not offered as a
   general real-data comparator.

The recorded numbers are single-run measurements on one machine; each
cell was measured once and variability across repeated runs was not
assessed. Recorded runs reused earlier outputs where noted: the pooled
real-AR(1) run reused the VOI/PPI inputs of an earlier exact run, and
the FAST-parent runs shared a precomputed real parent model. The
exact-only delivery (commit `f184cd66`) agrees with `main` to relative
max and L2 = 0 on all 22 fixture cases and the 3 PEB fixtures, but its
branch timings were not re-measured -- the speedup tables in
`REPORT.md` refer to the tested earlier versions.

### Synthetic single subject (optional; TMFC_simulations)

The `gppi/simulations/sim_sub01.csv` table was produced on subject 1 of the
paper simulations (designs 01/01v/12/13, 100 single-voxel ROIs,
`rng(1000*seed_id + subject)` with seed ids 1/201/12/13). The generation
follows the `TMFC_simulations` repository (`PIPELINE_01_BLOCK_DESIGN.m`
/ `PIPELINE_02_EVENT_RELATED_DESIGN.m`: SF=1, SNR=0.4, 0.2 s synaptic
delay, 3 dummy scans, TR=2 s, AR(1) parent). The generation and timing
harness is not part of this repository; the `.csv` and `.md` files are
the recorded evidence.
