# LSS evaluation

## Summary

`accelerated_lss` estimates the same LSS model as `main` faster using the
default `exact` method, including for rWLS parents. The optional `pooled`
method approximates the covariance calculation for AR(1), FAST and rWLS
parents: one sample covariance per run, per-trial ReML and whitening kept. `tmfc_LSS_after_FIR` uses a
batched solver and corrects repeated grand-mean scaling of FIR residuals.
The public interfaces and output files are unchanged. LSS runs serially:
parallel workers were slower in the measured configuration.

## Parameters and usage

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

## Results

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
[integration checks](sub-01/integration.md) and
[rWLS results](sub-01/rwls.md).

On five synthetic designs for one subject, exact max relative beta error per
trial was approximately 4.9e-15; pooled was approximately 5.9e-5 and BSC z
correlation rounded to 1.0000000. See [simulation results](simulations/sim_sub01.md).
After FIR, ten synthetic datasets were bitwise equal in float32. See
[after-FIR results](sub-01/after_fir.md).

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

### Ground truth on simulations

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

## Work in progress

FAST with `lss_reml_every > 1` needs further quality assessment: minimum
ROI series correlation was 0.862655 in a rare condition. The Python port is in
progress; whole-brain and timing comparisons remain open.

## How to test

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

## Appendix: methods and mathematics

### Per-trial model

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

### Noise model and ReML

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

### Exact estimation

Exact keeps the per-trial ReML and whitening of `main`. Shared across the
trials of a run: the loaded data, mask and global scaling, the high-pass
filter, the first-pass filtered data $`KY`$ (the first pass uses $`W = I`$)
and the data residuals on the constants/DCT basis. Only the order of
operations changes, so saved float32 betas agree (bitwise, or single-ULP
tails), and double-precision checks agree to 1.42e-12 absolute (AR(1)) or
4.56e-10 relative (end-to-end). Saved float32
agreement does not establish a double-precision bound; results depend on
platform and BLAS.

### Pooled estimation

Pooled replaces each trial's covariance by one per-run covariance
$`C_{y,0}`$ from a nuisance run model (all trials of interest replaced by
one regressor per run): $`C_{y,k} \approx C_{y,0}`$. For AR(1) and FAST
both the voxel pool $`\mathcal P_k`$ (the F selection at UFp) and the weights
are approximated; for rWLS the pool does not depend on the trial (UF = 0)
and only the weights are. Each trial still estimates its own $`h_k`$ and
$`W_k`$. Pooled changes the estimate; check its quality for your data.

### Warm start and interval

Warm start initializes ReML at the parent hyperparameters
$`h_\text{parent}`$ (for FAST, restricted to the parent's nonzero components
after model reduction). With `lss_reml_every = m`, ReML runs on every
$`m`$-th trial and intervening trials reuse the whitening of the nearest
estimated trial: $`W_k = W_{k'}`$. Both can change the result and carry no
acceptance gate.

### Rejected faster variants

Reusing the parent's fixed whitening matrix and fitting one ReML per run
were measured and not adopted. Fixed parent whitening produced
anticorrelated ROI series and sign flips on FAST (minimum ROI series correlation
−0.85, 16 flips): the parent whitening is estimated from a different
covariance pool and need not represent each trial. One ReML per run changed
rejected trial sets on FAST (BSC z correlation 0.992). Approximate rWLS
scan-weight variants also differed substantially on whole-brain data;
exact rWLS and run-pooled rWLS are the supported routes.

### After FIR and residual scaling

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
