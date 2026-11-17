# LSS after FIR: ds000008 sub-01

MATLAB R2025a, SPM25, serial, default BLAS threads; reference `main` with
the FIR residual scaling correction (branch `main_global_fixed`). Three runs, 546 scans, 387 trials and 84,821,499
in-mask voxel-trial entries. Each repetition ran in a separate MATLAB process
under `/usr/bin/time -l`.

| Method | Rep 1 s | Rep 2 s | Mean s | Speedup | Peak RSS rep 1 / 2 MiB | max rel L2 | max float32 ULP | Differing voxel-trial entries |
|---|---|---|---|---|---|---|---|---|
| corrected reference | 878.0 | 874.2 | 876.1 | 1 | 3860.7 / 3890.0 | reference | reference | reference |
| in_memory | 25.0 | 28.7 | 26.9 | 32.6 | 3860.0 / 3866.9 | 1.98e-13 | 1 | 6 |
| batched | 17.3 | 17.0 | 17.1 | 51.2 | 4678.8 / 4708.3 | 9.72e-11 | 1 | 7 |

Neither method had a mask mismatch. Synthetic tests on five designs and two
subjects (10 datasets, 1,760 trials) matched the corrected reference bitwise in
saved float32 betas. Synthetic and sub-01 comparisons use saved float32
outputs; they are not double-precision gates.

With `global='None'`, original `main` had run gSF 125.09, 125.82, 129.83.
Original betas divided by run gSF matched the corrected reference within 1 float32
ULP; BSC max |Δz| was 8.48e-9 for six voxel pairs in the first run and
condition. Real-data `Scaling`, complete ROI-to-ROI BSC and group analyses
were not run.
