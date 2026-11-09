#!/bin/bash
# Case Study 3: data for the read-based QKI splice maps (figure panel e).
# Run from a Mac terminal, in an environment with samtools:
#     conda activate rbp-maps
#     bash get_data.sh
# Everything is written under DATA_DIR (default: data/ in the case study folder):
#   encode_rmats/   ENCODE QKI knockdown rMATS tars (2)
#   liftover/       UCSC hg19ToHg38 chain
#   events/         GRCh38 event lists and the regions rbp-maps reads (prepare_events.py)
#   rbp_ark_eclip_skipper2/          per-region BAMs (+ library_sizes.tsv), ~1 GB
#   rbp_ark_eclip_skipper2/full/     the 7 full eCLIP BAMs of the two QKI datasets in the splice maps (3.7 GB)
# Safe to rerun: finished steps are skipped.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
DATA=${DATA_DIR:-$HERE/../data}
mkdir -p "$DATA/encode_rmats" "$DATA/liftover"
echo "Writing to $DATA"
for acc in ENCFF096LJD ENCFF103PLX; do
  [ -s "$DATA/encode_rmats/$acc.tar.gz" ] || curl -fL --retry 5 -o "$DATA/encode_rmats/$acc.tar.gz" \
    "https://www.encodeproject.org/files/$acc/@@download/$acc.tar.gz"
done
CHAIN="$DATA/liftover/hg19ToHg38.over.chain.gz"
[ -s "$CHAIN" ] || curl -fL --retry 5 -o "$CHAIN" \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
[ -s "$DATA/events/regions.bed" ] || python3 "$HERE/prepare_events.py" "$DATA/encode_rmats" "$CHAIN" "$DATA/events"
bash "$HERE/fetch_subset_bams.sh" "$HERE/manifest.tsv" "$DATA/events/regions.bed" "$DATA/rbp_ark_eclip_skipper2"
