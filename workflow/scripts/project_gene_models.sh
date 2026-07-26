#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  echo "Usage: $0 GENOME SOURCE_PROTEINS COMPETITIVE_FAI CONTIG_PREFIX OUTPUT_PREFIX EXTENSION THREADS OUTPUT_DIR" >&2
  exit 2
}
[[ $# -eq 8 ]] || usage

genome=$1
proteins=$2
fai=$3
contig_prefix=$4
output_prefix=$5
extension=$6
threads=$7
output_dir=$8
for file in "$genome" "$proteins" "$fai"; do
  [[ -s "$file" ]] || { echo "Missing input: $file" >&2; exit 1; }
done
mkdir -p "$output_dir"
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/gene_models.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT
miniprot -t "$threads" --gff "$genome" "$proteins" > "$tmpdir/projected.gff3"
python "$(dirname "$0")/build_quantseq_saf.py" \
  --gff "$tmpdir/projected.gff3" --fai "$fai" \
  --prefix "$contig_prefix" --already-prefixed 0 --extension "$extension" \
  --output-prefix "$output_prefix" --output-dir "$output_dir"
cp "$tmpdir/projected.gff3" "$output_dir/projected.gff3"
