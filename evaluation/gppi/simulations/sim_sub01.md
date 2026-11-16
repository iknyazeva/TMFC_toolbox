# Synthetic single subject (designs 01/01v/12/13, subject 1)

Settings: 100 single-voxel ROIs, gPPI analysis 2, strict mode; seeds
rng(1000*seed_id + 1): 01 -> 1001, 01v -> 201001, 12 -> 12001, 13 -> 13001.
Variants: exact (reference) and pooled; pooled warm start not run. The
design-12 exact reference is an earlier run of the same exact code on
inputs verified bit-identical to this subject's (non-paired reference);
its AUC row is recorded without a paired timing cell.

Subject-level AUC of the exact A-B contrast is near chance (0.49 - 0.54)
on one subject at SNR 0.4; agreement metrics are the primary result.
Agreement is computed on the A-B matrices (asymmetrical and symmetrical),
off-diagonal elements; flips count sign disagreements among |exact| > 0.1.

| design | AUC exact a/s | AUC pooled a/s | dAUC a/s | rel L2 a/s | Pearson r a/s | flips a/s | gPPI s exact/pooled |
|---|---|---|---|---|---|---|---|
| 01 | 0.5129 / 0.5127 | 0.5107 / 0.5105 | -0.0022 / -0.0022 | 0.0902 / 0.0757 | 0.9966 / 0.9977 | 1 / 0 | 191.3 / 105.6 |
| 01v | 0.4927 / 0.4909 | 0.4913 / 0.4892 | -0.0014 / -0.0017 | 0.1001 / 0.0791 | 0.9951 / 0.9970 | 4 / 0 | 197.3 / 107.7 |
| 12 | 0.5201 / 0.5202 | 0.5215 / 0.5228 | +0.0014 / +0.0026 | 0.1157 / 0.1145 | 0.9963 / 0.9974 | 2 / 0 | not paired / 90.9 |
| 13 | 0.5412 / 0.5451 | 0.5424 / 0.5467 | +0.0013 / +0.0016 | 0.1075 / 0.0991 | 0.9957 / 0.9970 | 3 / 2 | 250.9 / 127.3 |

a = asymmetrical, s = symmetrical. Machine-readable copy: sim_sub01.csv.
