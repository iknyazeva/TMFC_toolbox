# gPPI evaluation

## Summary

`accelerated_gppi` estimates the same gPPI and gPPI-FIR models as `main`
faster using the `exact` method (the only one in this delivery; a pooled
approximation existed pre-removal — see the historical sections). The acceleration covers
the whole path: VOI extraction, PPI regressor construction (fewer
`spm_PEB` deconvolution calls), and the per-ROI batch estimation in
`tmfc_gPPI` / `tmfc_gPPI_FIR` (batch template built once per subject, VOI
series read once, estimation loop restructured). Public interfaces and
output files are unchanged; runs are serial, without parallel workers.

On real data (ds000008 sub-01, 9 ROIs, serial) the exact method reproduces
`main` exactly as measured (relative max and L2 difference 0) while running
the gPPI branch 1.47x (AR(1)) and 1.20x (FAST) faster; on the 100-ROI
synthetic calibration the gPPI branch is 3.86x faster. Historically
(pre-removal), the pooled option was close to exact on AR(1) real data
(worst relative deviation 4.3e-3) with a further 1.77x (gPPI) and 1.88x
(gPPI-FIR) branch speedup, but deviated about 10% (relative L2) on the
single-voxel-ROI simulations and was pathological under FAST; see the
historical sections — pooled is NOT part of this delivery.

## Parameters and usage

This delivery is **exact-only** (owner decision 2026-11-17). The default
and only estimation method is `exact`: identical numerics to `main`,
faster. The former `pooled` and warm-start options were REMOVED from the
code; their full implementation is preserved in the branch history
(commit `3dc8823`, removed on top of it). Previously accepted option
values are rejected EXPLICITLY at the very top of both functions, before
any other configuration access, IO or SPM initialisation — no silent
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
covariance components a `None` parent does not have — pre-existing
analysis-2 limitation, unchanged).

## What the pooled method approximated (REMOVED — historical context)

The numbers in this section and in the pooled rows of the tables below
were measured on the PRE-REMOVAL code (commits `3dc8823`/`9d564bb`,
option still present). They are kept as research context and do NOT
validate the exact-only build of this working tree (which has not been
re-run).

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
analyzed runs — one pool and one `Cy` per subject, shared across seeds,
including cross-run blocks. The per-seed ReML estimation, whitening and
solving are kept (`tmfc_gPPI_pooled_estimate`), so only the covariance
input to ReML changes. The estimate is approximate by design.

## Timing (ds000008 sub-01, serial, one run per cell; seconds)

Exact rows were measured on pre-removal snapshots of the accelerated
sources (see AC1 audit); the exact numeric path is preserved verbatim in
this tree, but the exact-only build itself has not been re-timed.

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

Synthetic single-subject gPPI branch, seconds (see `simulations/`):

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

## Agreement

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
| Analysis 2, independent errors (`None`) | both methods | main's legacy estimation failure reproduced identically (accepted legacy restriction, owner-approved); NOT a numerical output agreement |
| Analysis 1 and 3 (seed-to-voxel) | both methods x AR(1) and `None` parents | numeric PASS, rel max/L2 = 0 |
| Analysis 1/3 x FAST or rWLS parents | NOT RUN (8 cells: 2 methods x 2 analyses x 2 noise models) | no coverage claim; shared estimation mechanics do not substitute for a run |

The executed matrix is narrower than the original acceptance criterion; no
full four-noise-model analysis-1/3 coverage is claimed.

Provenance of the exact-equivalence evidence: the direct `main`
comparisons above were executed on earlier validated source states (each
run under a per-start source snapshot). The final commit retains the
exact estimator blocks unchanged (static literal proof against the branch
base), and a replay with the final code reproduced the saved exact
outputs for the simulated design-01/13 and real AR(1) targets exactly
(relative L2 = 0). No fresh `main` comparison was run at the final commit
SHA: that part of the chain is code review plus snapshot transitivity,
not a same-SHA numerical run. The zero-error results above are the
measured values of the runs that were executed.

Pooled vs exact (report-only, no gate):

| Data | worst rel max | worst rel L2 | notes |
|---|---|---|---|
| Real AR(1), 4 matrix blocks | 4.30e-3 | 3.19e-3 | pooled; warm start 4.26e-3 |
| Synthetic Sub-1 (4 designs) | see `simulations/` | 0.076 - 0.116 (off-diagonal rel L2) | Pearson r >= 0.995 both forms; sign flips among |exact|>0.1: at most 4 |
| Real FAST | 8.6e+5 (gPPI), ~3e+2 (FIR) | same order | PATHOLOGICAL — pooled with FAST parents is guarded off |

## Synthetic single subject (designs 01/01v/12/13, subject 1)

See `simulations/sim_sub01.md`. Exact-method subject-level AUC of the A-B
contrast against the simulation ground truth is near chance (0.49 - 0.54)
on one subject at SNR 0.4, so the agreement metrics are the primary
result: pooled tracks exact with Pearson r >= 0.995 (both matrix forms,
all designs), off-diagonal relative L2 0.076 - 0.116, and at most 4 sign
flips among entries with |exact| > 0.1. AUC differences between pooled and
exact stay within 0.0026.

## Limitations

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
  pooled is not supported. Analysis 1/3 with rWLS (and FAST) parents were
  not run (see the coverage table).

## Why pooled agrees less on simulations

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
