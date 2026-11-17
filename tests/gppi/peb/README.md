# gPPI acceleration — PEB hoist: source + fixtures (EXECUTED: PASS)

Status: **EXECUTED 2026-11-16 — PASS (all 3 cases, all gates).** See
[`../../../evaluation/REPORT.md`](../../../evaluation/REPORT.md)
(Agreement) for the durable summary of this evidence. Summary: strict/repeated/centering_off → identical file names/counts
(10/12/8), finite, measured rel max and L2 = 0 (gate ≤1e-3),
spm_PEB calls 10→4 / 12→4 / 8→4; 2 MATLAB starts, 34 s cumulative.

## What changed (source diff)

`functions/tmfc_PEB_PPI.m` (+21/−6): the `spm_PEB` deconvolution
(`C = spm_PEB(Y,P); xn = xb*C{2}.E(1:N); xn = spm_detrend(xn)`) is hoisted
out of the (condition × ROI) loop into a per-ROI `xn_all` cell computed once
per run, plus a defensive `continue` for a run with an empty `jConds` list
(cannot happen with strict/repeated lists today). No approximation,
threshold, model, centering, convolution, timing, basis or output-name
changes. Invariants: `Y`, `xb`, `Hxb`, `X0`, `Q`, `W`, `P` depend only on
ROI × run, never on the condition.

## Files

| File | Purpose |
|---|---|
| `peb_ppi_fixture.m` | Small real-SPM parent (fmri_spec+fmri_est, AR(1), derivs [1 1], 2 runs, 12×12×8, 120 scans) + per-run VOI series; cases `strict`, `repeated`, `centering_off` |
| `run_peb_ppi_tests.m` | Runs reference (main `functions/`) vs candidate `tmfc_PEB_PPI`, gates + spm_PEB call counts via profiler |
| `README.md` | This file |

## Coverage matrix

| Case | session_mode | conds | bf | pmod | VOI style | centering |
|---|---|---|---|---|---|---|
| `strict` | strict | 3 + 2 across 2 runs, multi-cond sessions | 1, 2, 3 | B×Mod column | ROI_1 old (Y only), `ROI 2` Y+Yraw (Yraw=Y+1000, must be ignored; name has a space → regexprep) | with_mean_centering |
| `repeated` | repeated | 3 conds in every run | 1, 2 | B×Mod | same pair | with_mean_centering |
| `centering_off` | strict | 2 + 2, legacy cond_list (no `.pmod`/`.bf` fields → defaults) | 1 (default) | absent (catch path) | same pair | no_mean_centering |

Gates: identical file names/counts (== `expected_outputs`), `PPI.name`
equality, finiteness, relative max AND L2 ≤ 1e-3 (no ULP gate);
spm_PEB calls: reference = nROI × Σ(conds per run), candidate = nROI ×
runs-with-conds (strict 10→4, repeated 12→4, centering_off 8→4). Call
counting uses the MATLAB profiler `FunctionTable` (`FunctionName`,
`NumCalls`; the same method as the earlier phase-timing harness) — no path shadowing and no
`spm_PEB` stubs, so no recursion risk and no cleanup hazard.

## Reproduction (already executed; command for the record)

```matlab
matlab -batch "cd <this_checkout>/tests/gppi/peb; results = run_peb_ppi_tests( ...
    '/tmp/tmfc10_peb_out', ...
    '<main_toolbox>/functions', ...
    '<spm25>'); disp(results)"
```

Fresh `out_dir` outside the source tree; `ref_functions_dir` is the
read-only main toolbox checkout. Single case subset: pass e.g.
`{'strict'}` as 4th argument. No parfor involved.

## Known risks / non-goals

- Executed and passing on small synthetic fixtures (1 subject); not a
  cohort/ds000008 run; no ULP gate; agreement was measured as
  rel max and L2 = 0.
- LSS fixtures/tests are untouched; these tests live only in `tests/gppi/peb/`.
- The call-count check is structural (profiler-based); the numeric
  comparison against real SPM25 `spm_PEB` is the acceptance test and it
  passed on all cases.
