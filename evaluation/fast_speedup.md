# FAST acceleration measurements

Serial experiments used subject 3 for each synthetic design and sub-01 for
the real data. Baselines were pooled without warm start for design 02 and
sub-01, and exact for design 15 because the pooled baseline was unavailable.
The comparison for design 15 therefore differs from the other two.

`warm` initializes ReML with the parent parameters and FAST components;
`warm, interval 4` also reuses whitening between ReML estimates. These component
timings are distinct from the end-to-end results in the report.
Warm start is now available for pooled and is the default for FAST; interval 4
is experimental.

| Dataset | Option | Time s | Speedup | beta rel median / max | ROI corr median / min | BSC z corr min | max abs dz | sign flips | badRows Jaccard |
|---|---|---|---|---|---|---|---|---|---|
| 02 | warm | 139.6 | 3.60 | 0.0003 / 0.0007 | 1.00000 | 1.00000 | 0.0003 | 0 | 0.8–1 |
| 02 | warm, interval 4 | 42.8 | 11.74 | 0.0019 / 0.019 | 0.99999 | 0.99999 | 0.0028 | 0 | 1 |
| 15 | warm | 906.0 | 2.73 | 0.019 / 0.038 | 0.9997 | 0.99943 | 0.019 | 0 | 0.5–0.8 |
| 15 | warm, interval 4 | 248.0 | 9.96 | 0.019 / 0.038 | 0.9997 | 0.99943 | 0.019 | 0 | 0.5–0.8 |
| sub-01 vs pooled | warm | 28.3 | 1.88 | 0.0057 / 0.016 | 0.99998 / 0.982 | 0.9983 | 0.12 | 0 | 1 (–) |
| sub-01 vs pooled | warm, interval 4 | 18.6 | 2.86 | 0.0089 / 0.026 | 0.99993 / 0.875 | 0.9891 | 0.29 | 0 | 1 (–) |

For sub-01, minimum z correlation and max absolute z difference against
`main` were 0.9978 / 0.14 for warm and 0.9880 / 0.31 for warm with interval
4. The lowest correlation was in a rare run-2 `stop_failed` cell. The dash
marks a comparison where one rejected-trial set was empty. Synthetic ROI
correlations were recorded as one value, without a median/min pair.
Ground-truth AUC and group analyses were not measured.
