# Evaluation report

Combined report for the accelerated LSS and gPPI branches (one delivery line;
LSS history and gPPI history are both preserved below). All numbers, tolerances
and caveats are carried over verbatim from the per-method reports.

## Summary

### LSS

`accelerated_lss` estimates the same LSS model as `main` faster using the
default `exact` method, including for rWLS parents. The optional `pooled`
method approximates the covariance calculation for AR(1), FAST and rWLS
parents: one sample covariance per run, per-trial ReML and whitening kept. `tmfc_LSS_after_FIR` uses a
batched solver and corrects repeated grand-mean scaling of FIR residuals.
The public interfaces and output files are unchanged. LSS runs serially:
parallel workers were slower in the measured configuration.

### gPPI

`accelerated_gppi` estimates the same gPPI and gPPI-FIR models as `main`
faster using the `exact` method (the only one in this delivery; a pooled
approximation existed pre-removal -- see the historical sections). The
acceleration covers the per-ROI batch estimation in `tmfc_gPPI` /
`tmfc_gPPI_FIR` (batch initialisation and ROI-invariant design template
once per call/subject, VOI series read once, per-seed cloning) and the
`spm_PEB` deconvolution call pattern in `tmfc_PEB_PPI` (hoisted out of
the condition loop; call counts per fixture reduced from one per
condition-run to one per run). `tmfc_PPI` and `tmfc_VOI` themselves are
unchanged from the branch base (verified by diff); the PPI branch timing
difference below is consistent with the changed `spm_PEB` call pattern,
while the VOI branch differences may reflect run-to-run and cache
variability (inference from identical source; no repeated measurements
were taken; no speedup mechanism is claimed for VOI).
Public interfaces and output files are unchanged; runs are serial,
without parallel workers.

On real data (ds000008 sub-01, 9 ROIs, serial) the exact method reproduces
`main` exactly as measured (relative max and L2 difference 0) while running
the gPPI branch 1.47x (AR(1)) and 1.20x (FAST) faster; on the 100-ROI
synthetic calibration the gPPI branch is 3.86x faster. Historically
(pre-removal), the pooled option was close to exact on AR(1) real data
(worst relative deviation 4.3e-3) with a further 1.77x (gPPI) and 1.88x
(gPPI-FIR) branch speedup, but deviated about 10% (relative L2) on the
single-voxel-ROI simulations and was pathological under FAST; see the
historical sections -- pooled is NOT part of this delivery.

## LSS acceleration

| `tmfc.defaults` field | Values and behavior |
|---|---|
| `lss_method` | `'exact'` (default) or `'pooled'` (approximate). For rWLS parents, `'pooled'` shares one run covariance and keeps per-trial rWLS ReML and whitening. |
| `lss_reml_warm_start` | Logical; pooled only. Defaults to `true` for FAST parents and `false` otherwise; an explicit value overrides the default. Ignored for rWLS with a warning. |
| `lss_reml_every` | Positive integer; pooled only, default `1`. Values above `1` reuse whitening and are experimental. |
| `parallel` | `1` warns; LSS runs serially. |

Warm start and the interval are ignored for exact and rWLS estimation, with a
warning when explicitly supplied. Choose exact for agreement with `main`;
pooled trades exact agreement for speed. On FAST, inspect the warm-start
agreement for the conditions of interest; explicitly set it to `false` to
run pooled without warm start.

```matlab
tmfc.defaults.lss_method = 'exact';
tmfc_LSS(tmfc);

% Approximate run-pooled covariance; FAST warm start defaults to true.
tmfc.defaults.lss_method = 'pooled';
tmfc_LSS(tmfc);

% Experimental: reuse whitening between every fourth ReML estimation.
tmfc.defaults.lss_reml_every = 4;
tmfc_LSS(tmfc);
```

### Results

Whole-brain ds000008 sub-01, 380 trials, serial. Times are means of two
repetitions where two exist; the rWLS exact and pooled rows and the FAST
warm-start rows are single runs, so their speedups combine a single run
with a two-run mean.

| Noise model and method | Time, s | Speedup vs main | Agreement with main | Tolerance |
|---|---:|---:|---|---|
| None, main | 916.5 | 1× | reference | — |
| None, exact | 13.5 | 68× | max 1 float32 ULP | beta rel L2 ≤ 1e-3; names and masks identical |
| None, pooled | 13.9 | 66× | min ROI series corr 1.00000; min BSC z corr 1.00000 | min ROI series corr ≥ 0.999 |
| AR(1), main | 1362.3 | 1× | reference | — |
| AR(1), exact | 46.8 | 29× | beta rel L2 4.56e-10 (22 ULP) | beta rel L2 ≤ 1e-3; names and masks identical |
| AR(1), pooled; warm off | 22.4 | 61× | min ROI series corr 0.999996; min BSC z corr 1.00000 | min ROI series corr ≥ 0.999 |
| AR(1), pooled; warm on | 17.6 | 77× | min ROI series corr 0.999984; min BSC z corr 0.999992 | none — approximate, reported |
| AR(1), pooled; every 4 | 14.8 | 92× | min ROI series corr 0.999798; min BSC z corr 0.999967 | none — approximate, reported |
| FAST, main | 1561.4 | 1× | reference | — |
| FAST, exact | 78.7 | 20× | max 1 float32 ULP | beta rel L2 ≤ 1e-3; names and masks identical |
| FAST, pooled; warm off | 57.7 | 27× | min ROI series corr 0.99969; min BSC z corr 0.99997 | min ROI series corr ≥ 0.999 |
| FAST, pooled; warm on (default) | 35.1 | 44× | min ROI series corr 0.977497; min BSC z corr 0.997794 | none — approximate, reported |
| FAST, pooled; warm on, every 4 | 20.4 | 76× | min ROI series corr 0.862655; min BSC z corr 0.987975 | none — approximate, reported |
| rWLS, main | 1305.1 | 1× | reference | — |
| rWLS, exact | 67.1 | 19× | beta rel L2 4.40e-10 (462 ULP) | beta rel L2 ≤ 1e-3; names and masks identical |
| rWLS, pooled | 17.1 | 76× | min ROI series corr 0.999906; min BSC z corr 0.99996 | min ROI series corr ≥ 0.999 |
| After FIR, corrected reference | 876.1 | 1× | reference | — |
| After FIR, batched solver | 17.1 | 51× | beta rel L2 9.72e-11 (1 ULP) | beta rel L2 ≤ 1e-3; names and masks identical |

The AR(1) pooled warm-on (17.6 s) and every-4 (14.8 s) times are single
runs of a separate series; their speedups use the AR(1) main time of the
two-run timing series (1362.3 s). Against the pooled warm-off run of the
end-to-end checks (21.9 s) they are 1.24× and 1.48× faster.

In double precision, AR(1) exact agreed with `main` to a maximum absolute
beta difference of 1.42e-12; the tabulated relative value is from the
end-to-end check. The rWLS exact time is also end-to-end; an earlier
component-level run (87.3 s) was bitwise equal in saved float32, and data
reuse alone measured 89.0 s with zero difference. Warm start and interval
rows are approximation diagnostics without an acceptance gate; see
[integration checks](lss/sub-01/integration.md) and
[rWLS results](lss/sub-01/rwls.md).

On five synthetic designs for one subject, exact max relative beta error per
trial was approximately 4.9e-15; pooled was approximately 5.9e-5 and BSC z
correlation rounded to 1.0000000. See [simulation results](lss/simulations/sim_sub01.md).
After FIR, ten synthetic datasets were bitwise equal in float32. See
[after-FIR results](lss/sub-01/after_fir.md).

Beta rel is the trial L2 difference divided by reference L2; for whole-brain
approximate results the reported agreement uses two correlations with `main`:

- *ROI series corr*: for each condition and each of nine ROIs, the beta
  series (first eigenvariate of the ROI voxels, as in `tmfc_BSC`) is
  correlated across trials with the `main` series; trials rejected in
  either result are excluded from both. The table gives the minimum over
  ROIs and conditions.
- *BSC z corr*: correlation of the Fisher-z ROI-to-ROI matrices of the two
  results; minimum over conditions.

Both measure agreement with `main`, not accuracy against ground truth.
Measurements used a MacBook Pro M5 Pro, MATLAB R2025a, SPM25 and serial runs.

#### Ground truth on simulations

Paper simulations, 10 subjects, four event-related designs with 100 trials
per condition and TR 2 s. The ground truth is a binary ROI-by-ROI matrix of
edges whose connectivity differs between conditions A and B; only the
A−B contrast is evaluated. For each subject Δz = z_A − z_B (Fisher-z BSC);
a one-sample t-test across subjects gives the group statistic.

AUC is the probability that a random true edge receives a higher score than
a random null edge (the area under the ROC curve over all thresholds,
off-diagonal edges). The group score is |t|; the subject score is |Δz|.
0.5 is chance and 1 is perfect ranking. AUC needs no significance threshold;
with 10 subjects FDR-thresholded TPR is near zero, so AUC is the main
comparison.

| Design | AR(1) exact | AR(1) pooled | AR(1) one ReML per run | FAST pooled, warm on | After FIR |
|---|---:|---:|---:|---:|---:|
| ISI 6 s, variable HRF | 0.668 | 0.668 | 0.667 | 0.640 | 0.673 |
| ISI 6 s | 0.804 | 0.804 | 0.804 | 0.805 | 0.772 |
| ISI 8 s | 0.845 | 0.845 | 0.845 | 0.845 | 0.831 |
| ISI 12 s | 0.849 | 0.849 | 0.849 | 0.848 | 0.844 |

Group AUC of the A−B contrast, without trial rejection. AR(1) pooled and one
ReML per run match exact AR(1) to within 0.001 in group AUC, with and
without trial rejection. Pooled changed FDR significance (q = 0.05) of no
edge; one ReML per run changed at most one edge without rejection and up to
seven with rejection (ISI 12 s). One ReML per run is a research variant,
not a public option. FAST pooled with warm start matches AR(1) exact except
on the variable-HRF design (group AUC −0.028, −0.022 with rejection; mean
subject AUC −0.005, paired p ≈ 0.04 with rejection); exact FAST was not run, so
this difference cannot be split into noise-model and approximation effects.
After FIR compares LSS after FIR task regression with LSS without it on the
same images; it evaluates the analysis route, not the batched solver, whose
output equals the reference. With trial rejection, its mean subject AUC was
0.009–0.012 lower on the three designs with fixed HRF (paired p < 0.05).
Exact AR(1) matched `main` to a maximum absolute beta difference of
5e-13 on six subject–design pairs.

### Work in progress

FAST with `lss_reml_every > 1` needs further quality assessment: minimum
ROI series correlation was 0.862655 in a rare condition. The Python port is in
progress; whole-brain and timing comparisons remain open.

### How to test

Use `accelerated_lss` for the new methods and `main` for the original.
`main_global_fixed` is `main` with the FIR residual scaling correction.
Use separate writable project copies with the same parent models; LSS writes
beta files. [README](README.md) lists required inputs and commands.
Run `run_lss_tests(out_dir, ref_functions_dir, spm_dir)` in MATLAB;
`out_dir` must be disposable and `ref_functions_dir` a separate `main`
functions directory. Add `tests/lss` to the MATLAB path first.

The test suite uses its own pooled check, independent of ROI files: it
splits the brain into eight octants, averages the voxels of each octant
per trial, and correlates the upper triangles of the resulting 8×8
connectivity matrices of `main` and the tested method (*octant
connectivity corr*, gate ≥ 0.999). It is not the ROI series correlation
in the table above.

### Appendix: methods and mathematics

#### Per-trial model

For trial $`k`$, LSS fits $`y = X_k\beta + \varepsilon`$ with
$`X_k = [a_k, b_k, N]`$: the current-trial column, the pooled other-trials
column and the nuisance block (run constant or per-run constants plus the
parent confounds). With whitening $`W_k`$ and high-pass filtering applied to
data and design first,

```math
\widehat\beta_k = (X_k^\top W_k^\top W_k X_k)^{-1} X_k^\top W_k^\top W_k\, y,
```

solved through the rank-aware pseudo-inverse; only the first coefficient is
kept.

#### Noise model and ReML

$`\varepsilon \sim \mathcal N(0, \sigma^2 V)`$ with
$`V(h) = \sum_i h_i Q_i`$ for fixed covariance components $`Q_i`$. ReML
estimates $`h`$ from the weighted sample covariance of the filtered data over
the voxel pool $`\mathcal P_k`$ selected by the effects-of-interest F
contrast at the UFp threshold:

```math
C_{y,k} = \frac{1}{|\mathcal P_k|}\sum_{v \in \mathcal P_k}
\frac{\operatorname{tr}(RV)}{\operatorname{ResSS}_{k,v}}\; y_v y_v^\top,
\qquad W_k = V_k^{-1/2}.
```

AR(1) uses two components; FAST uses a larger dictionary with model
reduction; rWLS uses one component per scan (diagonal $`V`$, per-scan
weights) with a trial-independent pool (UF = 0). The design, the pool and
the residual sums depend on the trial, so $`h_k`$ and $`W_k`$ vary per trial.

#### Exact estimation

Exact keeps the per-trial ReML and whitening of `main`. Shared across the
trials of a run: the loaded data, mask and global scaling, the high-pass
filter, the first-pass filtered data $`KY`$ (the first pass uses $`W = I`$)
and the data residuals on the constants/DCT basis. Only the order of
operations changes, so saved float32 betas agree (bitwise, or single-ULP
tails), and double-precision checks agree to 1.42e-12 absolute (AR(1)) or
4.56e-10 relative (end-to-end). Saved float32
agreement does not establish a double-precision bound; results depend on
platform and BLAS.

#### Pooled estimation

Pooled replaces each trial's covariance by one per-run covariance
$`C_{y,0}`$ from a nuisance run model (all trials of interest replaced by
one regressor per run): $`C_{y,k} \approx C_{y,0}`$. For AR(1) and FAST
both the voxel pool $`\mathcal P_k`$ (the F selection at UFp) and the weights
are approximated; for rWLS the pool does not depend on the trial (UF = 0)
and only the weights are. Each trial still estimates its own $`h_k`$ and
$`W_k`$. Pooled changes the estimate; check its quality for your data.

#### Warm start and interval

Warm start initializes ReML at the parent hyperparameters
$`h_\text{parent}`$ (for FAST, restricted to the parent's nonzero components
after model reduction). With `lss_reml_every = m`, ReML runs on every
$`m`$-th trial and intervening trials reuse the whitening of the nearest
estimated trial: $`W_k = W_{k'}`$. Both can change the result and carry no
acceptance gate.

#### Rejected faster variants

Reusing the parent's fixed whitening matrix and fitting one ReML per run
were measured and not adopted. Fixed parent whitening produced
anticorrelated ROI series and sign flips on FAST (minimum ROI series correlation
−0.85, 16 flips): the parent whitening is estimated from a different
covariance pool and need not represent each trial. One ReML per run changed
rejected trial sets on FAST (BSC z correlation 0.992). Approximate rWLS
scan-weight variants also differed substantially on whole-brain data;
exact rWLS and run-pooled rWLS are the supported routes.

#### After FIR and residual scaling

After FIR the model sets `cvi = 'None'` and no high-pass filter: $`W = I`$.
The batched solver residualizes data and trial columns on the session
constant(s) $`N`$ once per run — $`r_N = (I - P_N)\,y`$ and
$`Z_k = (I - P_N)\,[a_k, b_k]`$ with $`P_N = N N^+`$ — and computes each trial
coefficient as the first row of $`\operatorname{pinv}(Z_k)\, r_N`$, using the
rank tolerance of `spm_sp` on column-normalized $`Z_k`$. All trials then need
one matrix product per data chunk. An in-memory per-trial alternative is
`tmfc_LSS_after_FIR_in_memory`. Against the corrected reference, saved
float32 synthetic outputs agreed bitwise and sub-01 differed by one ULP in
7 of 84,821,499 voxel-trial entries. Real-data `'Scaling'` mode and full
ROI-to-ROI BSC were not tested.

FIR residuals are already scaled by the parent model. The original routine
applies global scaling to these near-zero-mean residuals again. Under
`'None'` the spurious factor is constant within a run, so by linearity
$`\widehat\beta_\text{fixed} = \widehat\beta_\text{orig}/\mathrm{gSF}_\text{run}`$;
a common positive factor leaves within-run BSC correlations unchanged
(checked on six voxel pairs). Under `'Scaling'` the factor varies by scan
and distorts the time series and its correlations. Micro tests matched the
corrected result to unscaled OLS and the original result to OLS on
scan-scaled data.

## gPPI acceleration

### What is accelerated and how

(vs the branch base `3955863`)

Only three functions changed; `tmfc_PPI.m`, `tmfc_VOI.m` and the rest of
`functions/` are byte-identical to the base.

1. **Batch machinery built once instead of per-ROI** (`tmfc_gPPI.m`,
   `tmfc_gPPI_FIR.m`): `spm('defaults')` + `spm_get_defaults`
   registrations and `spm_jobman('initcfg')` run once per call
   (previously inside the per-ROI loop); the ROI-invariant part of the
   `fmri_spec` batch (timing, scans, conditions, nuisance regressors,
   HPF, basis set, masking/noise-model fields) is captured once per
   subject as a `base_spec` template and cloned per seed ROI, adding only
   the seed-dependent PPI/seed regressors and output directory.
   Batch semantics are preserved: SPM25 reads the registered defaults
   dynamically at run time, and nothing inside `initcfg`/`cfg_util`
   resets them (verified at the source level; the 22-case validation
   reproduces `main` exactly).
2. **VOI series read once per subject** (`tmfc_gPPI.m`, `tmfc_gPPI_FIR.m`):
   the per-run VOI `.mat` files are loaded once into `VOI_Yraw` (seed
   regressor, all analyses) and -- for the analysis modes that need the
   beta stack -- `Y` (filtered series), instead of one `load` per ROI in
   the estimation loop. Old-style VOI files (only `Y`, no `Yraw`) are
   still accepted through the documented fallback.
3. **PEB deconvolution hoisted out of the condition loop**
   (`tmfc_PEB_PPI.m`): the deconvolution inputs (`Y`, basis `xb`,
   neural `P`) depend only on ROI x run, not on the condition, so the
   per-ROI result is computed once and reused for every condition; a
   run with no selected conditions now skips the call entirely
   (defensive; cannot occur with the current strict/repeated lists).
   The estimator itself (`spm_PEB`) is unchanged -- the measured gain is
   the reduced call count (one per run instead of one per
   condition-run); no per-call speedup is claimed.

The noise-model branch (`spm_est_non_sphericity` vs the robust rWLS
helpers) and everything downstream of the beta computation (matrix
saving, beta-image copying, output names, batch `.mat` contents after
path normalisation) are unchanged; the exact code path is
preserved, and the agreement tables below show relative max and L2 = 0
on every verified case.

The removed `pooled` method reused the LSS machinery pattern
(`tmfc_LSS_pooled_cy`/`_estimate` conventions) for the covariance
computation; see the historical section.

### Parameters and usage

This delivery is **exact-only**. The default
and only estimation method is `exact`: identical numerics to `main`,
faster. The former `pooled` and warm-start options were REMOVED from the
code; their full implementation is preserved in the branch history
(commit `3dc8823`, removed on top of it). Previously accepted option
values are rejected EXPLICITLY at the very top of both functions, before
any other configuration access, IO or SPM initialisation -- no silent
fallback:

- `tmfc.defaults.gppi_method` present and not `'exact'` errors
  (`tmfc_gPPI:InvalidMethod` / `tmfc_gPPI_FIR:InvalidMethod`);
- `tmfc.defaults.gppi_reml_warm_start` present (any value) errors
  (`tmfc_gPPI:WarmStartRemoved` / `tmfc_gPPI_FIR:WarmStartRemoved`).

A basic call needs no option fields at all:

```matlab
tmfc_gPPI(tmfc);        % exact: same numerics as main, faster
tmfc_gPPI_FIR(tmfc);
```

Independent-error (`None`) parents still fail inside analysis 2 with the
same legacy error as `main` (SPM25 non-sphericity estimation requires
covariance components a `None` parent does not have -- pre-existing
analysis-2 limitation, unchanged).

### What the pooled method approximated (REMOVED -- historical context)

The numbers in this section and in the pooled rows of the tables below
were measured on the PRE-REMOVAL code (commits `3dc8823`/`9d564bb`,
option still present). They are kept as research context and refer to
the pre-removal code, not to the delivered exact-only build: that
build's accuracy was re-verified against `main` on the 22-case fixture
matrix and the PEB fixtures (relative max and L2 = 0), but its branch
timings were not re-measured.

In analysis 2 the exact path estimates the noise operator per seed ROI:
for each seed it strips the previous whitening, calls SPM's
`spm_est_non_sphericity` (ReML over the seed design's effects of interest
on the whole-brain pool), forms the whitening matrix and solves the
weighted least-squares problem. The pooled method replaces the per-seed
covariance ESTIMATE with one run/seed-pooled sample covariance per
subject: it builds a nuisance-only design (the seed design minus the PPI
and seed-regressor columns; FIR condition regressors stay), selects a
voxel pool with an F threshold, and accumulates `Y*Y'` of the selected,
variance-scaled raw time series over the concatenated timeline of all
analyzed runs -- one pool and one `Cy` per subject, shared across seeds,
including cross-run blocks. The per-seed ReML estimation, whitening and
solving are kept (`tmfc_gPPI_pooled_estimate`), so only the covariance
input to ReML changes. The estimate is approximate by design.

### Timing (ds000008 sub-01, serial, one run per cell; seconds)

Exact rows were measured on earlier validated source states of the
accelerated branch (per-start snapshots); the delivered exact-only build
(commit `f184cd66`) matches `main` to relative max and L2 = 0 on the 22-case validation
matrix, but its branch timings have not been re-measured -- treat the
speedups below as properties of the tested earlier versions.

Exact method vs `main`, per branch (no totals across branches):

| Branch | Parent | main, s | exact, s | ratio |
|---|---|---:|---:|---:|
| VOI | AR(1) | 10.42 | 5.62 | 1.86x |
| PPI | AR(1) | 1.51 | 1.15 | 1.31x |
| gPPI | AR(1) | 104.00 | 70.84 | 1.47x |
| gPPI-FIR | AR(1) | 105.37 | 68.89 | 1.53x |
| VOI | FAST | 9.5 | 4.5 | 2.1x |
| PPI | FAST | 1.6 | 1.3 | 1.2x |
| gPPI | FAST | 317.7 | 263.8 | 1.20x |
| gPPI-FIR | FAST | 388.7 | 334.1 | 1.16x |
| gPPI | AR(1), 100 single-voxel ROIs (synthetic calibration) | 592.0 | 153.3 | 3.86x |

Pooled vs exact (candidate), per branch, AR(1) real data:

| Branch | exact, s | pooled, s | pooled + warm start, s |
|---|---:|---:|---:|
| gPPI | 70.84 | 40.1 (1.77x) | 37.5 (1.89x) |
| gPPI-FIR | 68.89 | 36.7 (1.88x) | 43.1 (1.60x) |

Synthetic single-subject gPPI branch, seconds (see `gppi/simulations/`):

| Design | exact, s | pooled, s |
|---|---:|---:|
| 01 block | 191.3 | 105.6 |
| 01v variable HRF | 197.3 | 107.7 |
| 12 event ISI 2 s | not run paired (153.3 in the 100-ROI calibration) | 90.9 |
| 13 event ISI 4 s | 250.9 | 127.3 |

Not measured: real seed-to-voxel (analysis 1/3) timing (its accuracy is
fixture-covered; the timing is deferred), real-data timing for rWLS and
independent-error parents
(analysis-2 legacy limitation above), pooled on FAST real data beyond the
guarded run, pooled with rWLS, any multi-subject timing, parallel workers.

### Agreement

Exact vs `main` (gate: relative max AND L2 <= 1e-3; all measured 0):
14 of 14 executable fixture cases across noise models and analysis modes
(AR(1), FAST, rWLS, independent errors via analysis 1/3; analysis 1/3
image comparisons include geometry, NaN and Inf masks), the real-data
AR(1) and FAST comparisons above (VOI series, PPI regressors, asymmetrical
and symmetrical matrices for both methods), and the 100-ROI calibration
(observed rel max and L2 = 0 on all 9,900 finite matrix entries). The PEB deconvolution hoist
agrees with observed rel max and L2 = 0 on its fixtures (strict, repeated and centering-off) with
`spm_PEB` call counts reduced from 10/12/8 to 4 per fixture (2 ROIs x 2 runs; 2 calls per run).

Fixture coverage actually executed (numeric gate rel max AND L2 <= 1e-3,
all measured 0 unless stated):

| Fixture family | Coverage executed | Result |
|---|---|---|
| PEB deconvolution hoist | 3 cases (strict, repeated, centering-off) | rel max/L2 = 0; `spm_PEB` calls 10/12/8 -> 4 per fixture (2 ROIs x 2 runs) |
| Analysis 2 (ROI-to-ROI) | gPPI + gPPI-FIR x AR(1)/FAST/rWLS parents | numeric PASS, rel max/L2 = 0 |
| Analysis 2, independent errors (`None`) | both methods | main's legacy estimation failure reproduced identically (accepted legacy restriction); NOT a numerical output agreement |
| Analysis 1 and 3 (seed-to-voxel) | both methods x AR(1) and `None` parents | numeric PASS, rel max/L2 = 0 |
| Analysis 1/3 x FAST or rWLS parents | 8 cells (2 methods x 2 analyses x 2 noise models), executed in the final validation launch | numeric PASS, rel max/L2 = 0 |

Fixture coverage after the final validation launch spans all four noise
models in analysis 1/3 (AR(1), FAST, rWLS via fixtures plus independent
errors via the None a1/a3 cases) for both methods, in addition to the
analysis-2 matrix (independent-error parents excepted). The remaining scope limit is analysis 2 with
independent-error (`None`) parents: there the legacy estimation failure
of `main` is reproduced identically -- an accepted legacy restriction,
not a numerical output agreement.

Provenance of the exact-equivalence evidence: the historical real-data
and 100-ROI comparisons were executed on earlier validated source states
(each run under a per-start source snapshot); the performance numbers in
the timing tables refer to those source versions. The delivered exact-only
build (commit `f184cd66`) was validated directly against `main` in the
final validation launch: 22 fixture cases (14 regression + the 8
analysis-1/3 FAST/rWLS cells) plus the PEB fixtures and the option-guard
preflight, one serial process, rel max and L2 = 0 on every case, with
per-case source attribution (`which`) and pinned snapshots. Independent-error
parents in analysis 2 remain the accepted legacy failure (reproduced
identically, not a numerical output agreement). Real-data runs for
independent-error/rWLS parents and real seed-to-voxel remain deferred.
The zero-error results above are the measured values of the runs that
were executed.

Pooled vs exact (report-only, no gate):

| Data | worst rel max | worst rel L2 | notes |
|---|---|---|---|
| Real AR(1), 4 matrix blocks | 4.30e-3 | 3.19e-3 | pooled; warm start 4.26e-3 |
| Synthetic Sub-1 (4 designs) | see `gppi/simulations/` | 0.076 - 0.116 (off-diagonal rel L2) | Pearson r >= 0.995 both forms; sign flips among |exact|>0.1: at most 4 |
| Real FAST | 8.6e+5 (gPPI), ~3e+2 (FIR) | same order | PATHOLOGICAL -- pooled with FAST parents is guarded off |

### Synthetic single subject (designs 01/01v/12/13, subject 1)

See `gppi/simulations/sim_sub01.md`. Exact-method subject-level AUC of the A-B
contrast against the simulation ground truth is near chance (0.49 - 0.54)
on one subject at SNR 0.4, so the agreement metrics are the primary
result: pooled tracks exact with Pearson r >= 0.995 (both matrix forms,
all designs), off-diagonal relative L2 0.076 - 0.116, and at most 4 sign
flips among entries with |exact| > 0.1. AUC differences between pooled and
exact stay within 0.0026.

### Limitations

- Pooled deviates ~10% (relative L2) from exact on the synthetic
  single-voxel-ROI data vs 4.3e-3 on real data; the cause is the
  covariance input to ReML (see below), not voxel selection.
- Subject-level AUC is near chance on one subject; no multi-subject or
  group-level claims are made here.
- Warm start showed no measured accuracy or speed advantage on real AR(1)
  (4.26e-3 vs 4.30e-3; branch timings within noise) and is default-off.
- Pooled with FAST parents is numerically pathological at real scale
  (ReML objectives of order 1e9 during estimation; agreement rel max up to
  8.6e+5) and errored early with `PooledFastUnsupported` (pre-removal
  behaviour of the removed option; the current build rejects the option
  itself).
- The synthetic design-12 exact reference is a non-paired earlier run of
  the same exact code on identical inputs (verified bit-identical inputs),
  not a same-start pairing.
- rWLS parents: exact is validated on the analysis-2 fixtures above;
  pooled is not supported (removed). Analysis 1/3 with rWLS (and FAST)
  parents passed the fixture accuracy gate (8/8 cells, rel max and
  L2 = 0; see the coverage table); real-data timing for these parents
  is not measured.

### Why pooled agrees less on simulations

Diagnostic on three targets: simulated designs 01 and 13 (10 seeds each)
and real sub-01 AR(1) (all 9 seeds). Per seed it rebuilt the exact and
pooled voxel pools and `Cy` from one data load, re-estimated both methods,
and ran a counterfactual: the pooled estimation fed the seed's own exact
`Cy`. The deviation
metric is the rel L2 over the PPI rows of the per-seed beta matrices
(distinct from the final A-B matrix agreement reported for the synthetic
subject above). Validity first: the re-estimated betas reproduce the saved exact
matrices exactly and the pooled ones to <= 2.4e-14.

| Target | voxel pools (exact = pooled) | dev_pooled rel L2, PPI rows, per seed (mean / max) | counterfactual rel L2, PPI rows (mean / max) | `Cy` rel Frobenius difference |
|---|---|---|---|---|
| Sim design 01 | identical, 100 voxels | 0.08374 / 0.21326 | 0 / 0 (exact, all seeds) | 1.0 |
| Sim design 13 | identical, 100 voxels | 0.12689 / 0.29647 | 0 / 0 (exact, all seeds) | 1.0 |
| Real sub-01 | identical, 219,177 voxels | 0.00197 / 0.00331 | 0.00053 / 0.00061 | 0.013 - 0.024 |

- Voxel selection is not the cause: the pools are identical everywhere
  (on simulations every in-mask voxel is selected by both F-tests).
- The covariance input carries the whole gap on the simulations and most
  of it on the real data: with the seed's own `Cy` the pooled estimation
  reproduces exact betas exactly on simulations (counterfactual rel L2 0)
  and cuts the real-data rel-L2 gap 3.75x (mean) / 5.47x (max), with a
  nonzero real residual (mean 5.3e-4, max 6.1e-4) left unresolved below. Across seeds the beta gap tracks the
  whitening distance (Spearman 0.98 / 0.85 / 0.83, n = 10 / 10 / 9).
- Why the regimes differ: `Cy` weights each voxel by `sqrt(trRV/ResSS)`
  of the residualising design, which is nuisance-only for pooled and the
  full seed design (seed and PPI columns) for exact. On real data the two
  residualisations nearly coincide (`Cy` difference 1-2%). On simulations
  the difference equals the full `Cy` norm to printed precision for every
  seed, and the exact ReML hyperparameters are much larger than pooled,
  with some pooled fits blowing up to ~1e150. This is consistent with,
  but does not directly measure, the following mechanism: each ROI is one
  voxel, the seed's own voxel is fitted almost perfectly by its seed
  regressor, and its near-zero `ResSS` gives it a weight that dominates
  the exact `Cy`. The exact noise model on these simulations is then
  effectively estimated from the seed voxel alone, while pooled uses all
  100 voxels. The ~10% gap is therefore a property of single-voxel-ROI
  data and is not expected at real ROI sizes. Per-voxel weights and `Cy`
  norms were not recorded, so this mechanism is not verified.
- Cross-run concatenation is not implicated: the 3-run real data and the
  1-run simulations show the same pattern.
- Residual real-data gap after the counterfactual (~6e-4): not resolved
  further.
