# ds000008 sub-01: noise model comparisons

Serial, 380 trials. Two timings are repetitions in reversed order. Relative
beta error is absolute difference divided by voxel beta-series SD. ROI and
BSC correlations are minima over condition/run cells; badRows compares
rejected-trial sets. Pooled FAST used warm start off. Fixed parent whitening,
one ReML per run and approximate rWLS weights are research comparisons.

| Noise model | Method | rel median / p99 / max | ROI corr min | BSC z corr min | flips | badRows Jaccard | time rep 1 / 2, s | speedup |
|---|---|---|---|---|---|---|---|---|
| none | main | reference | — | — | — | — | 907.0 / 926.0 | 1 |
| none | in-memory reference | 0 / 0 / 0 | 1 | 1 | 0 | 1 | 30.5 / 30.7 | 30 |
| none | exact | 1ulp / 1ulp / 1ulp | 1 | 1 | 0 | 1 | 13.4 / 13.5 | 68 |
| none | pooled warm off | 0 / 0 / 1.3e-7 | 1.00000 | 1.00000 | 0 | 1.00 | 13.6 / 14.1 | 66 |
| none | fixed parent whitening | 0 / 0 / 2.9e-7 | 1.00000 | 1.00000 | 0 | 1.00 | 15.2 / 15.4 | 60 |
| none | one ReML per run | 0 / 0 / 2.9e-7 | 1.00000 | 1.00000 | 0 | 1.00 | 15.2 / 15.5 | 60 |
| fast | main | reference | — | — | — | — | 1492.2 / 1630.6 | 1 |
| fast | in-memory reference | 0 / 0 / 0 | 1 | 1 | 0 | 1 | 305.5 / 314.1 | 5 |
| fast | exact | 1ulp / 1ulp / 1ulp | 1 | 1 | 0 | 1 | 76.6 / 80.8 | 20 |
| fast | pooled warm off | 4.3e-4 / 3.8e-3 / 1.8e-2 | 0.99969 | 0.99997 | 0 | 1.00 | 53.2 / 62.2 | 27 |
| fast | fixed parent whitening | 4.4e-2 / 1.9e-1 / 6.8e-1 | -0.851 | 0.350 | 16 | 0.67 | 17.3 / 17.2 | 90 |
| fast | one ReML per run | 1.2e-2 / 5.4e-2 / 1.8e-1 | 0.92308 | 0.99248 | 0 | 0.67 | 18.2 / 18.4 | 85 |
| rwls | main | reference | — | — | — | — | 1229.8 / 1380.5 | 1 |
| rwls | exact | 0 / 0 / 0 | 1 | 1 | 0 | 1 | 85.5 / 89.1 | 15 |
| rwls | run nuisance weights | 7.7e-2 / 4.9e-1 / 2.4 | 0.832 | 0.94368 | 0 | 0.50 | 15.2 / 16.8 | 82 |
| rwls | parent weights | 7.4e-2 / 4.7e-1 / 2.3 | 0.325 | 0.80296 | 16 | 0.50 | 14.5 / 16.2 | 85 |
| rwls | run all-trial weights | 1.3e-1 / 7.4e-1 / 2.8 | -0.922 | 0.695 | 24 | 0.50 | 15.5 / 17.8 | 78 |

With independent errors W is the identity, so the methods coincide mathematically; written float32 files can differ by an ULP.
