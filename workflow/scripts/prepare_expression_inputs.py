#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
import yaml


def read_fasta(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    descriptions: dict[str, str] = {}
    sequences: dict[str, str] = {}
    current = ""
    chunks: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    sequences[current] = "".join(chunks)
                header = line[1:].strip()
                parts = header.split(maxsplit=1)
                current = parts[0]
                descriptions[current] = parts[1] if len(parts) > 1 else ""
                chunks = []
            elif current:
                chunks.append(line.strip())
    if current:
        sequences[current] = "".join(chunks)
    return descriptions, sequences


def clean_description(text: str) -> str:
    return re.sub(r"\s*\[[^\]]+\]\s*$", "", text).strip()


def source_accession(target: object) -> str:
    text = str(target or "").strip()
    return text.split()[0] if text and text.lower() != "nan" else ""


def write_fasta(records: list[tuple[str, str, str]], path: Path) -> None:
    with path.open("wt", encoding="utf-8", newline="\n") as handle:
        for identifier, description, sequence in records:
            if not sequence:
                continue
            handle.write(f">{identifier} {description}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare fungal expression matrices and functional gene sets.")
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--sample-qc", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--usable-genes", type=Path, required=True)
    parser.add_argument("--gene-metadata", type=Path, required=True)
    parser.add_argument("--source-proteins", type=Path, required=True)
    parser.add_argument("--functional-terms", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-count", type=int, default=5)
    parser.add_argument("--min-samples", type=int, default=3)
    args = parser.parse_args()

    counts = pd.read_csv(args.counts, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str, keep_default_na=False)
    sample_qc = pd.read_csv(args.sample_qc, sep="\t")
    mapping = pd.read_csv(args.mapping, sep="\t")
    background = pd.read_csv(args.background, sep="\t", dtype=str, keep_default_na=False)
    usable = pd.read_csv(args.usable_genes, sep="\t")["gene_id"].astype(str).tolist()
    metadata = pd.read_csv(args.gene_metadata, sep="\t", dtype=str, keep_default_na=False)
    descriptions, sequences = read_fasta(args.source_proteins)
    terms = yaml.safe_load(args.functional_terms.read_text(encoding="utf-8"))
    compiled = {name: [re.compile(pattern, re.I) for pattern in patterns] for name, patterns in terms.items()}

    fus = manifest[manifest["pathogen"].str.lower().eq("fus")].copy()
    all_samples = fus["sample_id"].tolist()
    qc_map = sample_qc.set_index("sample_id")["pass"].astype(str).str.lower().isin(["true", "1", "yes"])
    primary_samples = [sample for sample in all_samples if bool(qc_map.get(sample, False))]
    if len(all_samples) != 24 or len(primary_samples) < 18:
        raise RuntimeError(f"Unexpected inoculated sample counts: all={len(all_samples)}, primary={len(primary_samples)}")

    background_set = set(background.loc[background["exclude"].astype(str).str.lower().isin(["true", "1", "yes"]), "gene_id"])
    inclusive_genes = [gene for gene in counts.index if (counts.loc[gene, all_samples] >= args.min_count).sum() >= args.min_samples]
    primary_genes = [gene for gene in usable if gene in counts.index]
    if len(primary_genes) < 100:
        raise RuntimeError(f"Too few primary genes: {len(primary_genes)}")

    annotation = metadata.copy()
    if "gene_id" not in annotation:
        raise RuntimeError("Gene metadata lacks gene_id")

    def text_column(name: str) -> pd.Series:
        if name in annotation.columns:
            return annotation[name].fillna("").astype(str)
        return pd.Series("", index=annotation.index, dtype=str)

    annotation["source_protein"] = text_column("Target").map(source_accession)
    annotation["product"] = annotation["source_protein"].map(lambda key: clean_description(descriptions.get(key, "")))
    annotation["description_text"] = (
        text_column("product") + " " + text_column("description") + " " + text_column("Name")
    ).str.strip()

    module_rows: list[dict[str, str]] = []
    for row in annotation.itertuples(index=False):
        text = str(getattr(row, "description_text", ""))
        for module, patterns in compiled.items():
            matches = [pattern.pattern for pattern in patterns if pattern.search(text)]
            if matches:
                module_rows.append({"module": module, "gene_id": row.gene_id, "matched_terms": ";".join(matches)})
    modules = pd.DataFrame(module_rows)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    counts.loc[inclusive_genes, all_samples].to_csv(args.output_dir / "counts_all.tsv", sep="\t", index=True, index_label="gene_id", lineterminator="\n")
    counts.loc[primary_genes, all_samples].to_csv(args.output_dir / "counts_primary.tsv", sep="\t", index=True, index_label="gene_id", lineterminator="\n")
    mapping_columns = [column for column in ["sample_id", "target_unique_rpm", "post_fastp_reads"] if column in mapping.columns]
    fus = fus.merge(mapping[mapping_columns], on="sample_id", how="left", validate="one_to_one")
    if "target_unique_rpm" not in fus.columns or fus["target_unique_rpm"].isna().any():
        raise RuntimeError("Mapping summary lacks target_unique_rpm for one or more inoculated samples")
    fus["ft_unique_rpm_input"] = fus["target_unique_rpm"]
    fus.to_csv(args.output_dir / "samples_all.tsv", sep="\t", index=False, lineterminator="\n")
    fus[fus["sample_id"].isin(primary_samples)].to_csv(args.output_dir / "samples_primary.tsv", sep="\t", index=False, lineterminator="\n")
    annotation.to_csv(args.output_dir / "gene_annotation.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")
    modules.to_csv(args.output_dir / "functional_gene_sets.tsv", sep="\t", index=False, lineterminator="\n")

    records = []
    annotation_index = annotation.set_index("gene_id")
    for gene in primary_genes:
        source = str(annotation_index.loc[gene, "source_protein"]) if gene in annotation_index.index else ""
        records.append((gene, str(annotation_index.loc[gene, "product"]) if gene in annotation_index.index else "", sequences.get(source, "")))
    write_fasta(records, args.output_dir / "projected_source_proteins.faa")


if __name__ == "__main__":
    main()
