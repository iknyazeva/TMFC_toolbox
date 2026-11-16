# Integrated sub-01 checks

Public functions on `accelerated_lss` were compared with saved references;
file names and masks matched. The acceptance criterion was maximum trial
relative beta L2 error ≤ 1e-3. Maximum ULP is reported separately.

| Route / reference | Maps | Time, s | max relative beta L2 | max ULP |
|---|---:|---:|---:|---:|
| exact / main | 380 | 61.746 | 4.56e-10 | 22 |
| pooled / saved pooled reference | 380 | 21.889 | 0 | 0 |
| after-FIR batched / corrected reference | 387 | 18.366 | 9.72e-11 | 1 |

The AR(1) pooled route had minimum mean-ROI beta-series correlation
0.999999687 against main. FAST warm start and warm start with interval 4
matched the saved experimental outputs for 380 maps (zero relative error and
zero ULP). Approximate results against main have no quality acceptance gate.
