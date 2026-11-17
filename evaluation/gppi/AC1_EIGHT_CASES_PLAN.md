# AC1 eight-case completion + exact-only validation — ONE launch packet (NOT RUN)

Prepared 2026-11-17 in the clean worktree `TMFC-10-gppi-clean`
(`clean/accelerated_gppi`). **Execution state = the owner's removal
commit** (the exact-only edits in this tree, committed by the owner after
review; commits `3dc8823` + `9d564bb` preserve the full pooled/warm
history locally — only until/unless history is rewritten, which nobody
plans). The launch refuses to run until that commit exists and matches
`EXPECTED_PIN`. Reference main: `TMFC-17-ihb-main` @ `3fb347b`, read-only.
**Zero MATLAB starts were performed; nothing below is claimed PASS.**

## 1. What this single launch validates

(a) The exact-only build after the pooled/warm removal: `tmfc_gPPI.m` /
`tmfc_gPPI_FIR.m` keep the exact arms verbatim (indent-insensitive diff
vs `9d564bb`: only pooled-code deletions + the option guards); the five
`tmfc_gPPI_pooled_*.m` helpers and `tests/gppi/options/` are removed (no
other code referenced them; LSS pooled helpers untouched).
(b) The eight missing AC1 fixture cells: `fast_s2v`, `fast_a1`,
`rwls_s2v`, `rwls_a1`, `fir_fast_s2v`, `fir_fast_a1`, `fir_rwls_s2v`,
`fir_rwls_a1` (real FAST parents; real rWLS parents via the clean-dir
`parent_wls` path; seeds 12/13 per the s2v/a1 convention).
(c) Post-removal regression of EVERY previously executed fixture mode
(all 14 old cases): `ar1`, `fir_ar1`, `fast`, `fir_fast`, `rwls`,
`fir_rwls` (a2 core), `ar1_a1`, `fir_ar1_a1` (a1), `ar1_s2v`,
`fir_ar1_s2v` (a3), `none_s2v`, `none_a1`, `fir_none_s2v`,
`fir_none_a1` (None a1/a3). Total: **22 exact-matrix cases**.
(d) The PEB hoist (`tests/gppi/peb`, unchanged source) and the guard
preflight `test_gppi_option_guards.m` (throws on failure; verifies the
`InvalidMethod`/`WarmStartRemoved` rejections — the guards sit at the top
of both functions, before any tmfc access, IO or SPM init, so firing
there proves rejection precedes all heavy work; the defaults-free path is
validated end-to-end by the fixture cases, not by this preflight).

One MATLAB process, serial, stop-on-first-error (runner rethrow), no
retry without fresh approval. Gates unchanged: rel max AND rel L2 <= 1e-3
per case, batch equality after root tokenization, identical
counts/names/shapes/NaN patterns, NaN-only-diagonal matrices, beta-image
geometry + NaN/Inf masks + non-empty finite mask (a1/a3).

## 2. The exact command (executable as-is; no hidden edits; no rm -rf)

```bash
set -euo pipefail
CLEAN=/Users/irina/Work/Projects/active/TMFC/worktrees/TMFC-10-gppi-clean
MAIN=/Users/irina/Work/Projects/active/TMFC/worktrees/TMFC-17-ihb-main
OUT=/Users/irina/Work/Projects/active/TMFC/work/tmfc10_gppi/exact_matrix/M7_exact_only_validation
WATCH=/Users/irina/Work/Projects/active/TMFC/work/tmfc10_gppi/exact_matrix/exact_matrix_watch.py
export SPM_DIR=/Users/irina/Documents/MATLAB/spm
: "${EXPECTED_PIN:?set EXPECTED_PIN to the owner-approved removal commit SHA}"

# (0) pin validation: HEAD must BE the approved commit and the tree clean
[ "$(git -C "$CLEAN" rev-parse HEAD)" = "$EXPECTED_PIN" ] || { echo "HEAD != EXPECTED_PIN"; exit 5; }
[ -z "$(git -C "$CLEAN" status --porcelain)" ] || { echo "worktree not clean - refuse"; exit 5; }

# (1) fail if any output already exists (NO rm -rf anywhere in this packet)
[ ! -e "$OUT" ] || { echo "output exists: $OUT"; exit 6; }
mkdir -p "$OUT/run" "$OUT/logs"

# (2) snapshot EXACTLY the existing sources (candidate tests+functions incl.
#     robust rWLS helpers; main comparison targets; the watchdog itself)
( cd "$CLEAN" && shasum -a 256 tests/gppi/exact/*.m tests/gppi/peb/*.m \
    functions/tmfc_gPPI.m functions/tmfc_gPPI_FIR.m functions/tmfc_PEB_PPI.m \
    functions/denoising/tmfc_spm_rwls_spm.m \
    functions/denoising/tmfc_spm_rwls_reml.m \
    functions/denoising/tmfc_spm_rwls_est_non_sphericity.m ; \
  cd "$MAIN" && shasum -a 256 functions/tmfc_gPPI.m functions/tmfc_gPPI_FIR.m ; \
  shasum -a 256 "$WATCH" ) | sort -k2 > "$OUT/run/sources.sha256"

# (3) preflight pgrep -fl MATLAB (user's hard requirement): STOP on any
#     remaining line after excluding ONLY matlab-mcp-core-server. pgrep
#     exits 1 on no matches -> '|| true' keeps set -e safe; the decision
#     uses the filtered log, and is rechecked immediately before launch.
pgrep -fl MATLAB > "$OUT/run/pgrep_before.txt" 2>&1 || true
BAD="$(grep -v 'matlab-mcp-core-server' "$OUT/run/pgrep_before.txt" || true)"
[ -z "$BAD" ] || { echo "REFUSING: other MATLAB-related processes:"; printf '%s\n' "$BAD"; exit 4; }

# (4) owned watchdog launch (usage per its source: <label> <timeout_s>
#     <cleanup_reserve_s> <command...>; owned PGID, TERM -> grace ->
#     per-PID SIGKILL re-verified by pgid; never kills foreign processes)
export CALIB_LEDGER_DIR="$OUT/run" CALIB_SOURCES_SHA="$OUT/run/sources.sha256" \
       CALIB_TIME_FILE="$OUT/logs/M7.time" CALIB_PREFLIGHT_EVIDENCE="$OUT/run/pgrep_before.txt"
printf '{"unit":"AC1 8 cases + exact-only validation","pin":"%s","cases":22}\n' "$EXPECTED_PIN" > "$OUT/logs/ctx.json"
export CALIB_CONTEXT="$OUT/logs/ctx.json"
BAD2="$(pgrep -fl MATLAB 2>/dev/null | grep -v 'matlab-mcp-core-server' || true)"
[ -z "$BAD2" ] || { echo "REFUSING at T-0:"; printf '%s\n' "$BAD2"; exit 4; }
python3 "$WATCH" M7_exact_only_validation 1300 150 \
  /usr/bin/time -l -o "$OUT/logs/M7.time" \
  /Applications/MATLAB_R2025a.app/bin/matlab -batch "\
    addpath('$CLEAN/tests/gppi/exact');\
    test_gppi_option_guards('$CLEAN/functions');\
    r1 = run_gppi_exact_tests('$OUT/out','$MAIN/functions','$SPM_DIR',{\
      'ar1','fir_ar1','fast','fir_fast','rwls','fir_rwls',\
      'ar1_a1','fir_ar1_a1','ar1_s2v','fir_ar1_s2v',\
      'none_s2v','none_a1','fir_none_s2v','fir_none_a1',\
      'fast_s2v','fast_a1','rwls_s2v','rwls_a1',\
      'fir_fast_s2v','fir_fast_a1','fir_rwls_s2v','fir_rwls_a1'});\
    addpath('$CLEAN/tests/gppi/peb');\
    r2 = run_peb_ppi_tests('$OUT/peb','$MAIN/functions','$SPM_DIR');\
    disp('M7_DONE')" \
  > "$OUT/logs/M7.log" 2>&1
```

The guard preflight throws inside MATLAB on any unexpected outcome (its
return value is not a boolean to ignore). Read results from the SAVED
files, not the console: `$OUT/out/results.json` (22 per-case records with
`which` source attribution per variant), `$OUT/peb/results.mat`, ledger
`$OUT/run/ledger_M7_exact_only_validation.json`. Pass = every case's rel
gates, batch_equal and structural checks true; PEB rel max/L2 <= 1e-3 (the recorded value was 0) with reduced
`spm_PEB` counts (4 per fixture: 2 ROI x 2 runs).

## 3. Ceiling (actual ledger anchors; startup counted ONCE; RSS = estimate)

Observed watchdog walls: `out_final1` 44.1 s (ar1), `out_cont04` 46.1 s
(fir_ar1), M1 `fir_fast` 34.4 s, M5 `rwls+fir_rwls` 61.0 s, M2
`ar1_a1+fir_ar1_a1` 69.0 s, M3 `ar1_s2v+fir_ar1_s2v` 71.0 s, M6 four
None a1/a3 160.6 s, PEB unit ~20 s — each of these seven ledgers embeds
one ~20-25 s MATLAB startup; in a SINGLE launch that startup is paid
once. Marginal (startup-subtracted) sum for the 14 regression cases +
PEB + guard preflight ≈ (44+46+34+61+69+71+160+20) − 7×22 + 5 ≈ **370 s
observed**. The eight NEW cells never ran — EXTRAPOLATION from
same-analysis (M2/M3 ≈ 24 s/case marginal) and same-noise (M5, rWLS
parents ≈ 9 s/case extra) anchors: FAST a1/a3 ≈ 25-30 s each, rWLS a1/a3
≈ 35-45 s each ≈ **250-300 s extrapolated**. Single-launch total ≈ 370 +
275 + 25 ≈ **670 s**; proposed external watchdog ceiling **1300 s**
(≈1.9x — honest headroom for the rWLS a1/a3 per-seed robust path, the
main uncertainty). Cleanup reserve 150 s per watchdog semantics.
**RSS**: the watchdog logs peak RSS but enforces no limit — expect
<= 2 GB (estimate from M1-M6/PEB peaks 0.9-1.61 GB), stated as estimate.
One MATLAB launch total; on any failure STOP (ledger + evidence
preserved); continuation needs fresh approval — no retry, no source edits.

## 4. Static verification performed in this preparation (nothing executed)

Indent-insensitive diff of both functions vs `9d564bb`: only deletions +
the corrected top-of-function guards (properly quoted strings; guard now
precedes the first nargin handling — verified by line inspection in both
files). if/end balance signatures identical to HEAD (bracket-aware
scanner, `1:end` indexing handled). ASCII on all edited files;
`git diff --check` clean; 22-case list cross-checked against both
fixture/meta switches; guard test contains no mkdir/load/SPM calls and
throws instead of returning a boolean. References to removed
options/helpers updated in README.md (gPPI section exact-only),
evaluation/gppi/README.md and REPORT.md (pooled = historical context
only). Historical evidence and LSS sections untouched.

## 5. Blockers / approvals needed

Owner: review + create the removal commit, then communicate its SHA as
`EXPECTED_PIN` (the launch verifies HEAD == EXPECTED_PIN and a clean
tree; the sources.sha256 snapshot is evidence, and HEAD validation is the
commit check). Irina: ceiling approval for the single 1300 s (+150 s
cleanup) MATLAB launch. Publication / Done — separate approvals. Nothing
in this packet runs MATLAB.
