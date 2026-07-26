#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  echo "Usage: $0 MANIFEST BAM_DIR GENE_MODEL_DIR MODEL_PREFIX EXTENSION THREADS OUTPUT_DIR" >&2
  exit 2
}
[[ $# -eq 7 ]] || usage
manifest=$1
bam_dir=$2
models=$3
model_prefix=$4
extension=$5
threads=$6
output_dir=$7
mkdir -p "$output_dir"
mapfile -t samples < <(awk -F'\t' 'NR>1 {print $1}' "$manifest")
(( ${#samples[@]} > 0 )) || { echo "No samples in $manifest" >&2; exit 1; }
bams=()
for sample in "${samples[@]}"; do
  bam="$bam_dir/${sample}.bam"
  [[ -s "$bam" ]] || { echo "Missing BAM: $bam" >&2; exit 1; }
  samtools quickcheck -q "$bam"
  bams+=("$bam")
done
for annotation in native "ext${extension}"; do
  saf="$models/${model_prefix}_${annotation}.saf"
  [[ -s "$saf" ]] || { echo "Missing SAF: $saf" >&2; exit 1; }
  for strand in 0 1 2; do
    output="$output_dir/counts_${annotation}_s${strand}.txt"
    featureCounts -T "$threads" -F SAF -a "$saf" -o "$output" -s "$strand" --primary "${bams[@]}"
  done
done
