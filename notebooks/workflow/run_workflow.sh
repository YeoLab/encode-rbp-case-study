#!/bin/bash
# Run the Case Study 3 workflow notebooks 01-07 in order, saving their outputs in place.
# Start it from this folder with the encode-rbp-case-study environment active:  bash run_workflow.sh
# FROM=05 bash run_workflow.sh starts at notebook 05. Logs go to <data>/run_logs/workflow_*.log.
set -uo pipefail
ulimit -n 10240 2>/dev/null || ulimit -n 4096 2>/dev/null || true   # rbp-maps opens many files; macOS allows 256
cd "$(dirname "$0")"
DATA=${CASE_STUDY_DATA:-../../data}
LOGS="$DATA/run_logs"; mkdir -p "$LOGS"
command -v jupyter >/dev/null || { echo "jupyter not found: activate the encode-rbp-case-study environment"; exit 1; }
KEEP_AWAKE=""; command -v caffeinate >/dev/null && KEEP_AWAKE="caffeinate -i"
for nb in 01_gene_search 02_experiment_detail 03_download 04_shrna_rnaseq 05_differential_splicing 06_dpsi_per_exon 07_integrate; do
  [[ "${nb:0:2}" < "${FROM:-01}" ]] && continue
  echo "$(date '+%H:%M') running $nb"
  if ! $KEEP_AWAKE jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 \
       --ExecutePreprocessor.kernel_name=python3 "$nb.ipynb" > "$LOGS/workflow_$nb.log" 2>&1; then
    echo "$nb failed: see $LOGS/workflow_$nb.log and the notebook's last output"; exit 1
  fi
done
echo "$(date '+%H:%M') done: ../../figures/workflow/CS3_PICALM_QKI_workflow.png"
