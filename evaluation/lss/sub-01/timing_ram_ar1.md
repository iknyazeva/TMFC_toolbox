# ds000008 sub-01: AR(1), serial time and RAM

MacBook Pro M5 Pro, MATLAB R2025a and SPM25; 380 trials. The table reports
two-run mean time. RAM is peak MATLAB process-tree memory and process RSS,
respectively. Pooled used warm start off. Whitening reuse and one ReML per
run are research comparisons. Worker-process RAM was not measured.

| Method | Time s | s/trial | Speedup | Peak process-tree RAM GB | Peak RSS GB |
|---|---|---|---|---|---|
| main | 1362.3 | 3.59 | 1 | 5.5-5.8 | 4.80 |
| in-memory reference | 111.1 | 0.29 | 12.3 | 4.4-5.1 | 3.24-3.99 |
| exact | 46.8 | 0.12 | 29.1 | 4.9-5.2 | 3.78-4.08 |
| pooled | 22.4 | 0.06 | 60.8 | 5.5 | 4.36-4.37 |
| whitening reuse every 4 | 24.1 | 0.06 | 56.7 | 4.3 | 4.31-4.35 |
| one ReML per run | 16.9 | 0.04 | 80.6 | 4.1 | 4.07-4.08 |

With 18 workers, exact took 618.1 s and pooled 427.4 s; workers used one computational thread, while the serial process used the default BLAS threads.
