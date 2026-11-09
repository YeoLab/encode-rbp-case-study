#!/bin/bash
# Run Case Study 3 notebooks 01-05 one after another, saving their outputs in place.
# Start it from this folder in a terminal where the encode-rbp-case-study environment is active:
#     bash run_all.sh
# It first waits for notebook 00 (open in Jupyter Lab) to finish, i.e. for ../figures/data/picalm_locus.json.
# Logs: <data>/run_logs/. Stops at the first notebook that fails. Keeps the Mac awake while running.
# Optional: RBPMAPS_SIGTEST=zscore bash run_all.sh   (quick maps; default is the permutation test, ~5 h)
#           FROM=02 bash run_all.sh                    (start at notebook 02, skipping finished ones)
#           FROM=00 bash run_all.sh                    (run notebook 00 here too instead of in Jupyter Lab)
set -uo pipefail
ulimit -n 10240 2>/dev/null || ulimit -n 4096 2>/dev/null || true   # macOS defaults to 256 open files; rbp-maps opens many
cd "$(dirname "$0")"
DATA=${CASE_STUDY_DATA:-../data}
LOGS="$DATA/run_logs"; mkdir -p "$LOGS"
STATUS="$LOGS/status.tsv"
# log MESSAGE: prints MESSAGE with a UTC timestamp and appends it to <data>/run_logs/status.tsv
log() { printf "%s\t%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$STATUS"; }
command -v jupyter >/dev/null || { echo "jupyter not found: activate the encode-rbp-case-study environment first"; exit 1; }
for t in rmats.py rmats2sashimiplot samtools bedtools bedGraphToBigWig; do command -v $t >/dev/null || { echo "$t not found: wrong environment?"; exit 1; }; done

log "start (sigtest=${RBPMAPS_SIGTEST:-permutation}, source=${SE_SOURCE:-rerun}, from=${FROM:-01})"
while [ "${FROM:-01}" != "00" ] && [ ! -s ../figures/data/picalm_locus.json ]; do
  log "waiting for notebook 00 to finish (figures/data/picalm_locus.json)"; sleep 120
done
[ "${FROM:-01}" = "00" ] || log "notebook 00 finished"
FROM=${FROM:-01}
for nb in 00_screen_and_locus 01_download_shrna_bams 02_run_rmats 03_panel_c_sashimi 04_eclip_rbp_maps 05_window_maps_optional; do
  [[ "${nb:0:2}" < "$FROM" ]] && { log "skipping $nb (FROM=$FROM)"; continue; }
  log "running $nb"
  if caffeinate -i jupyter nbconvert --to notebook --execute --inplace \
       --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=python3 "$nb.ipynb" > "$LOGS/$nb.log" 2>&1; then
    log "finished $nb"
  else
    log "FAILED $nb (see $LOGS/$nb.log and the notebook's last output)"; exit 1
  fi
done
log "all done"
