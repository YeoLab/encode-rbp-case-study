#!/bin/bash
# Download eCLIP BAMs from RBP-ARK one at a time and make the small per-region BAMs that rbp-maps reads.
#
#   bash fetch_subset_bams.sh <manifest.tsv> <regions.bed> <out_dir>
#
# manifest.tsv: dataset <tab> path on RBP-ARK's file server <tab> bytes <tab> sha256
# For each BAM: download (resumable) to <out_dir>/full/ -> check SHA-256 -> keep read 2 if paired -> keep reads
# overlapping regions.bed -> <out_dir>/<name>.regions.bam (+ .bai). The full library size (samtools view -c -F 260,
# the count rbp-maps uses for RPM) goes to <out_dir>/library_sizes.tsv so maps can be normalized as if from the
# full BAM. Full BAMs are kept (KEEP_FULL=1, default); KEEP_FULL=0 deletes each after subsetting.
# Needs curl, shasum and samtools. Rerunning skips finished files.
set -euo pipefail
MANIFEST=$1; REGIONS=$2; OUT=$3
HOST=${RBPARK_FILES:-https://g-f8db23.0ed28.75bc.data.globus.org}
KEEP_FULL=${KEEP_FULL:-1}
command -v samtools >/dev/null || { echo "samtools not found: activate a conda env that has it (e.g. conda activate rbp-maps)"; exit 1; }
mkdir -p "$OUT/full"
SIZES="$OUT/library_sizes.tsv"
[ -s "$SIZES" ] || printf "bam\tdataset\tpaired\tmapped_F260\tkept_reads\tsha256_ok\n" > "$SIZES"

# free_bytes: prints the free space, in bytes, on the drive that holds OUT
free_bytes() { df -k "$OUT" | awk 'NR==2 {print $4 * 1024}'; }

n=0; total=$(grep -c . "$MANIFEST")
while IFS=$'\t' read -r ds path bytes sha <&3; do
  [ -n "$ds" ] || continue
  n=$((n + 1))
  f=$(basename "$path"); full="$OUT/full/$f"; sub="$OUT/${f%.bam}.regions.bam"
  if [ -s "$sub.bai" ] && grep -q "^$f	" "$SIZES"; then echo "[$n/$total] $f: done already"; continue; fi
  have=0; [ -f "$full" ] && have=$(wc -c < "$full" | tr -d ' ')
  if [ "$have" != "$bytes" ]; then
    need=$((bytes + 400000000 - have))
    if [ "$(free_bytes)" -lt "$need" ]; then
      echo "Not enough free disk for $f: need $((need / 1000000)) MB, have $(( $(free_bytes) / 1000000 )) MB. Free space and rerun."; exit 1
    fi
    echo "[$n/$total] $f ($((bytes / 1000000)) MB)"
    curl -fL -# --retry 5 --retry-delay 30 -C - -o "$full" "$HOST$path"
  else
    echo "[$n/$total] $f: already downloaded"
  fi
  got=$(shasum -a 256 "$full" | cut -d' ' -f1)
  if [ "$got" != "$sha" ]; then
    mv "$full" "$full.bad"; echo "  SHA-256 mismatch ($got); kept as $f.bad. Delete it and rerun."; exit 1
  fi
  # Paired BAMs (ENCORE/ENCORE2) keep read 2 only, as rbp-maps' walkthrough does; ENCODE3 genome_R2 BAMs are R2 already
  paired=$( (samtools view -f 1 "$full" 2>/dev/null || true) | head -1 | wc -l | tr -d ' ')
  FLT=""; [ "$paired" = "1" ] && FLT="-f 128"
  mapped=$(samtools view -c -F 260 $FLT "$full")
  samtools view -b $FLT -L "$REGIONS" -o "$sub" "$full"
  samtools index "$sub"
  kept=$(samtools view -c -F 260 "$sub")
  printf "%s\t%s\t%s\t%s\t%s\tyes\n" "$f" "$ds" "$paired" "$mapped" "$kept" >> "$SIZES"
  if [ "$KEEP_FULL" = "1" ]; then [ -s "$full.bai" ] || samtools index "$full"; else rm -f "$full"; fi
  echo "  paired=$paired mapped=$mapped kept=$kept ($(du -h "$sub" | cut -f1))"
done 3< "$MANIFEST"
echo "All $total BAMs done: $OUT"
