# Synthetic designs: subject 1 versus `main`

MATLAB R2025a, SPM25, serial; AR(1) parent, UFp 1 and HPF 128 s. Values
summarize each design across conditions and raw/clean BSC. Relative beta error
is the maximum trial L2 difference divided by reference L2. Beta-series
correlation is the minimum over simulated ROI series. BSC uses Fisher-z;
flips count reference |z| > 0.1, and badRows Jaccard compares rejected trials.
ΔAUC is the method minus main A–B raw AUC. Elapsed time and maximum elementwise
|Δz| were not recorded. Values are rounded as shown.

| Design | Trials | Scans | Method | max abs beta diff | max relative L2 | beta corr min | BSC z corr min | z mean bias max | flips | badRows Jaccard | ΔAUC |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 02 | 200 | 706 | exact | 4.386e-13 | 4.372e-15 | 1.0000000 | 1.0000000 | 7.54e-16 | 0 | 1.00 | +0.0e+00 |
| 02 | 200 | 706 | pooled | 2.197e-03 | 2.090e-05 | 1.0000000 | 1.0000000 | 2.69e-06 | 0 | 1.00 | +1.5e-06 |
| 02v | 200 | 706 | exact | 3.784e-13 | 3.753e-15 | 1.0000000 | 1.0000000 | 2.73e-16 | 0 | 1.00 | +0.0e+00 |
| 02v | 200 | 706 | pooled | 3.177e-03 | 2.202e-05 | 1.0000000 | 1.0000000 | 8.36e-06 | 0 | 1.00 | -4.6e-06 |
| 12 | 200 | 308 | exact | 3.659e-13 | 3.806e-15 | 1.0000000 | 1.0000000 | 1.41e-16 | 0 | 1.00 | +0.0e+00 |
| 12 | 200 | 308 | pooled | 4.260e-03 | 3.470e-05 | 1.0000000 | 1.0000000 | 2.65e-06 | 0 | 1.00 | -6.4e-06 |
| 15 | 200 | 1317 | exact | 4.743e-13 | 4.902e-15 | 1.0000000 | 1.0000000 | 4.52e-16 | 0 | 1.00 | +0.0e+00 |
| 15 | 200 | 1317 | pooled | 1.630e-03 | 1.465e-05 | 1.0000000 | 1.0000000 | 1.91e-06 | 0 | 1.00 | +0.0e+00 |
| 17 | 80 | 287 | exact | 3.570e-13 | 3.513e-15 | 1.0000000 | 1.0000000 | 8.59e-16 | 0 | 1.00 | +0.0e+00 |
| 17 | 80 | 287 | pooled | 6.446e-03 | 5.939e-05 | 1.0000000 | 1.0000000 | 1.12e-05 | 0 | 1.00 | -1.3e-06 |
