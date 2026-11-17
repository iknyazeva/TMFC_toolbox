# Synthetic timing: subject 3 calibration

Serial, default BLAS threads, 200 trials. Fixed parent whitening and one ReML per run are research comparisons. Subject-1 timings were not recorded.

| design (scans) | exact AR | pooled AR | fixed parent whitening AR | one ReML per run AR | fixed parent whitening AR UFp 0.001 | exact FAST | pooled FAST | fixed parent whitening FAST | one ReML per run FAST |
|---|---|---|---|---|---|---|---|---|---|
| 02 (706)   | 57.7 | 40.8 | 13.8 | 14.4 | 16.5 | 558.6 | 501.8 | 17.1 | 19.8 |
| 02v (706)  | 35.5 | 35.7 | 11.5 | 11.6 | 11.4 | 552.8 | 548.7 | 18.0 | 20.0 |
| 14 (~1000) | 56.0 | 55.5 | 14.4 | 14.7 | 14.1 | 861.0 | 853.3 | 20.9 | 25.9 |
| 15 (~1400) | 119.2 | 119.2 | 20.8 | 21.8 | 21.9 | 2470.0 | not measured* | not measured | not measured |

Seconds. The guard stopped the remaining FAST runs for design 15; those entries were not measured.
