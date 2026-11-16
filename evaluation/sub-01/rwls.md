# ds000008 sub-01: rWLS LSS

Source: end-to-end measurements on `accelerated_lss`,
single serial runs, MATLAB R2025a, SPM25, default BLAS threads.
Reference: the rWLS branch of `main` (1229.8 / 1380.5 s over two runs; see
[modes_none_fast_rwls](modes_none_fast_rwls.md)). 380 trials, whole brain.
Acceptance: max beta rel L2 <= 1e-3 for exact; min ROI corr >= 0.999 for
pooled.

| Method | Time, s | Agreement with main | Notes |
|---|---:|---|---|
| exact | 67.1 | max beta rel L2 4.40e-10; max 462 ULP | reuses the filtered data across trials |
| exact, data reuse only | 89.0 | zero difference; 0 ULP | isolates the data-reuse saving |
| exact, component benchmark | 87.3 | saved float32 bitwise equal | component protocol, two runs |
| pooled | 17.1 | max beta rel L2 8.12e-4; min ROI corr 0.999906 | one run covariance; per-trial rWLS ReML and W |

Exact kept the same relative error and ULP count against the earlier
bitwise component run; file names and masks matched both references. The previous
sub-01 gate of at most 1 ULP does not hold for the end-to-end exact route;
ULP counts are reported, not gated. Warm start and the interval are ignored
for rWLS with a warning.

Pooled diagnostics (nine conditions, worst cell; beta error relative to the
voxelwise beta-series SD): 8.97e-4 median / 2.14e-2 p99 / 1.17e-1 max; BSC
z correlation min 0.99996052; max abs dz 0.01048; 0 sign flips;
rejected-trial Jaccard 1 in the seven conditions where the metric is
defined. The min ROI corr >= 0.999 criterion passed.

Micro fixtures (single, outliers, concat): exact zero difference and zero
ULP against `main`; pooled outliers/concat minimum ROI correlation
0.9999598 / 0.9999972, BSC z correlation 0.999995 / 0.999999, no sign
flips, rejected-trial Jaccard 1. Degenerate zero-residual cases and
multi-chunk robustness beyond concat are not fully qualified. Single
subject and parent; no ground-truth or group-level assessment.
