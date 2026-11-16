# ds000008 sub-01 timing and agreement (serial, single runs)

Exact method (this branch) vs main, per branch; gate for exact agreement
is relative max AND L2 <= 1e-3 on identical inputs.

## timing_ar1.csv / timing_fast.csv

Branch-level walls, AR(1) and FAST parents (9 ROIs), and the 100-ROI
synthetic calibration (AR(1), design 12, subject 1).

## agreement_exact.md

Real-data AR(1) and FAST: exact vs main — identical file names, counts,
NaN masks; relative max and L2 = 0 on VOI series (54), PPI regressors
(81), gPPI and gPPI-FIR asymmetrical/symmetrical matrices (9 + 9 each).

Approximations on AR(1) (report-only, no gate): pooled worst rel max
4.30e-3 / rel L2 3.19e-3; pooled + warm start 4.26e-3 / 3.19e-3; VOI and
PPI inputs shared bit-identically with the exact run.

Pooled on FAST is guarded off: in the recorded run before the guard its
agreement was pathological (gPPI rel max 8.6e+5, FIR ~3e+2) and pooled
gPPI offered no speed advantage; the option errors early with
`PooledFastUnsupported`.
