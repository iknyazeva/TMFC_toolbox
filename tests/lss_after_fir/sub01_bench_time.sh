#!/bin/bash
# Run each benchmark in its own MATLAB process. GPL-3.0-or-later.
# Required environment: MATLAB_BIN, TMFC_BENCH_OUT, TMFC_PARENT_SPM, SPM_DIR.
set -eu
: "${MATLAB_BIN:?Set MATLAB_BIN to the MATLAB executable}"
: "${TMFC_BENCH_OUT:?Set TMFC_BENCH_OUT to the output directory}"
: "${TMFC_PARENT_SPM:?Set TMFC_PARENT_SPM to the parent SPM.mat}"
: "${SPM_DIR:?Set SPM_DIR to the SPM directory}"
TESTS="$(cd "$(dirname "$0")" && pwd)"
export TMFC_BENCH_OUT TMFC_PARENT_SPM SPM_DIR
mkdir -p "$TMFC_BENCH_OUT/time_logs"
matlab_batch() {
    "$MATLAB_BIN" -batch "addpath('$TESTS'); run_lss_after_fir_sub01_bench($1,getenv('TMFC_BENCH_OUT'),getenv('TMFC_PARENT_SPM'),getenv('SPM_DIR'))"
}
matlab_batch "'fir',[],[]"
for m in main in_memory batched; do
    for r in 1 2; do
        log="$TMFC_BENCH_OUT/time_logs/${m}_${r}.time"
        /usr/bin/time -l -o "$log" "$MATLAB_BIN" -batch "addpath('$TESTS'); run_lss_after_fir_sub01_bench('run','$m',$r,getenv('TMFC_BENCH_OUT'),getenv('TMFC_PARENT_SPM'),getenv('SPM_DIR'))"
    done
done
matlab_batch "'scaling_fix',[],[]"
matlab_batch "'report',[],[]"
