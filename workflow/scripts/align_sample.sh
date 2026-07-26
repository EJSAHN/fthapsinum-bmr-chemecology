#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  echo "Usage: $0 FASTQ SAMPLE INDEX CONTIG_MAP TARGET_BED TARGET_GROUP MIN_MAPQ THREADS OUT_BAM OUT_SPECIES OUT_SUMMARY OUT_FASTP_JSON OUT_COVERAGE" >&2
  exit 2
}
[[ $# -eq 13 ]] || usage

fastq=$1
sample=$2
index=$3
contig_map=$4
target_bed=$5
target_group=$6
min_mapq=$7
threads=$8
out_bam=$9
out_species=${10}
out_summary=${11}
out_fastp_json=${12}
out_coverage=${13}

for file in "$fastq" "$contig_map" "$target_bed"; do
  [[ -s "$file" ]] || { echo "Missing input: $file" >&2; exit 1; }
done
[[ -s "${index}.1.ht2" || -s "${index}.1.ht2l" ]] || { echo "Missing HISAT2 index: $index" >&2; exit 1; }

mkdir -p "$(dirname "$out_bam")" "$(dirname "$out_species")" "$(dirname "$out_summary")" "$(dirname "$out_fastp_json")" "$(dirname "$out_coverage")"
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/${sample}.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT
trimmed="$tmpdir/${sample}.fastq.gz"
full_bam="$tmpdir/${sample}.competitive.bam"
fastp_html="$tmpdir/${sample}.fastp.html"

fastp \
  --in1 "$fastq" --out1 "$trimmed" \
  --qualified_quality_phred 20 --unqualified_percent_limit 40 \
  --length_required 30 --trim_poly_x --thread "$threads" \
  --json "$out_fastp_json" --html "$fastp_html"

hisat2 -x "$index" -U "$trimmed" \
  --rna-strandness F --dta --max-intronlen 50000 -k 5 \
  --new-summary --summary-file "$out_summary" -p "$threads" \
  | samtools sort -@ "$threads" -m 768M -o "$full_bam" -

samtools quickcheck -q "$full_bam"
python "$(dirname "$0")/count_species.py" \
  --bam "$full_bam" --contig-map "$contig_map" --sample "$sample" \
  --min-mapq "$min_mapq" --output "$out_species"

prefix="$tmpdir/${sample}.target"
mosdepth -n --thresholds 1,3,5 --by "$target_bed" "$prefix" "$full_bam"
cp "${prefix}.thresholds.bed.gz" "$out_coverage"

mapfile -t contigs < <(awk -F'\t' -v group="$target_group" 'NR>1 && $2==group {print $1}' "$contig_map")
(( ${#contigs[@]} > 0 )) || { echo "No contigs for group $target_group" >&2; exit 1; }
samtools view -@ "$threads" -b -q "$min_mapq" -F 2308 -o "$out_bam" "$full_bam" "${contigs[@]}"
samtools quickcheck -q "$out_bam"
samtools index -@ "$threads" "$out_bam"
