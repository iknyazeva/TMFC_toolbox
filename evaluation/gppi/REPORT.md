# gPPI evaluation

## Summary

`accelerated_gppi` estimates the same gPPI and gPPI-FIR models as `main`
faster using the default `exact` method, and adds an optional `pooled`
approximation for the ROI-to-ROI noise estimation. The acceleration covers
the whole path: VOI extraction, PPI regressor construction (fewer
`spm_PEB` deconvolution calls), and the per-ROI batch estimation in
`tmfc_gPPI` / `tmfc_gPPI_FIR` (batch template built once per subject, VOI
series read once, estimation loop restructured). Public interfaces and
output files are unchanged; runs are serial, without parallel workers.

On real data (ds000008 sub-01, 9 ROIs, serial) the exact method reproduces
`main` bit-for-bit (relative max and L2 difference 0) while running the
gPPI branch 1.47x (AR(1)) and 1.20x (FAST) faster; on the 100-ROI
synthetic calibration the gPPI branch is 3.86x faster. The optional pooled
method is close to exact on AR(1) real data (worst relative deviation
4.3e-3) with a further 1.77x (gPPI) and 1.88x (gPPI-FIR) branch speedup,
but deviates about 10% (relative L2) on the single-voxel-ROI simulations;
it is guarded off for FAST parents, where its agreement was pathological.

## Parameters and usage

| `tmfc.defaults` field | Values and behavior |
|---|---|
| `gppi_method` | `'exact'` (default) or `'pooled'`. Shared by `tmfc_gPPI` and `tmfc_gPPI_FIR`. Invalid values error. |
| `gppi_reml_warm_start` | Logical, default `false`; pooled only. Warm-started ReML initialisation. Setting it with `exact` warns and ignores it. |

Guards (errors with specific identifiers, checked before any estimation):
`pooled` is supported for ROI-to-ROI analysis (analysis 2) with AR(1)
parents only — pooled with analysis 1/3 errors (`PooledUnsupported`);
pooled with a FAST parent errors (`PooledFastUnsupported`; numerically
pathological at real scale, see below); pooled with an rWLS parent errors
(`PooledRWLSUnsupported`); pooled with concatenated sessions errors
(`PooledConcatUnsupported`); pooled with independent errors (`None`)
fails inside analysis 2 with the same legacy error as `main` and `exact`
(SPM25 non-sphericity estimation requires covariance components that a
`None` parent does not have — a pre-existing limitation of the analysis-2
path, unchanged by this branch).

```matlab
tmfc.defaults.gppi_method = 'exact';      % default: same numerics as main, faster
tmfc_gPPI(tmfc);

tmfc.defaults.gppi_method = 'pooled';     % approximation, AR(1) + analysis 2 only
tmfc.defaults.gppi_reml_warm_start = false;
tmfc_gPPI(tmfc);
```

## What the pooled method approximates

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

Not measured: real-data timing for rWLS and independent-error parents
(analysis-2 legacy limitation above), pooled on FAST real data beyond the
guarded run, pooled with rWLS, any multi-subject timing, parallel workers.

## Agreement

Exact vs `main` (gate: relative max AND L2 <= 1e-3; all measured 0):
14 of 14 executable fixture cases across noise models and analysis modes
(AR(1), FAST, rWLS, independent errors via analysis 1/3; analysis 1/3
image comparisons include geometry, NaN and Inf masks), the real-data
AR(1) and FAST comparisons above (VOI series, PPI regressors, asymmetrical
and symmetrical matrices for both methods), and the 100-ROI calibration
(bitwise on all 9,900 finite matrix entries). The PEB deconvolution hoist
agrees bitwise on its fixtures (strict, repeated and centering-off) with
`spm_PEB` call counts reduced from 10/12/8 to 4 per run.

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
  8.6e+5) and errors early with `PooledFastUnsupported`.
- The synthetic design-12 exact reference is a non-paired earlier run of
  the same exact code on identical inputs (verified bit-identical inputs),
  not a same-start pairing.
- rWLS parents: exact is validated (fixtures); pooled is not supported.

## Why pooled agrees less on simulations

Diagnostic on three targets: simulated designs 01 and 13 (10 seeds each)
and real sub-01 AR(1) (all 9 seeds). Per seed it rebuilt the exact and
pooled voxel pools and `Cy` from one data load, re-estimated both methods,
and ran a counterfactual: the pooled estimation fed the seed's own exact
`Cy`. Validity first: the re-estimated betas reproduce the saved exact
matrices exactly and the pooled ones to <= 2.4e-14.

| Target | voxel pools (exact = pooled) | rel L2 pooled vs exact, PPI rows (mean / max) | counterfactual (exact `Cy`) | `Cy` rel Frobenius difference |
|---|---|---|---|---|
| Sim design 01 | identical, 100 voxels | 0.11 / 0.26 | 0 (all seeds) | 1.0 |
| Sim design 13 | identical, 100 voxels | 0.13 / 0.27 | 0 (all seeds) | 1.0 |
| Real sub-01 | identical, 219,177 voxels | 0.0019 / 0.0030 | 0.00059 / 0.00072 | 0.013 - 0.024 |

- Voxel selection is not the cause: the pools are identical everywhere
  (on simulations every in-mask voxel is selected by both F-tests).
- The covariance input carries the whole gap: with the seed's own `Cy`
  the pooled estimation reproduces exact betas exactly on simulations and
  cuts the real-data gap about 3.3x. Across seeds the beta gap tracks the
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
