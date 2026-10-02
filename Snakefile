configfile: "config.yaml"

import re
from pathlib import Path

import pandas as pd

shell.executable("/bin/bash")

SAMPLES = pd.read_csv(config["samples"], sep="\t", dtype=str, keep_default_na=False)
SETS = pd.read_csv(config["analysis_sets"], sep="\t", dtype=str, keep_default_na=False)
REFERENCES = pd.read_csv(config["references"], sep="\t", dtype=str, keep_default_na=False)

SAMPLE_INFO = SAMPLES.set_index("sample_id").to_dict("index")
REFERENCE_INFO = REFERENCES.set_index("label").to_dict("index")


def members(set_name):
    return SETS.loc[SETS["set_name"].eq(set_name), "sample_id"].tolist()


FUSARIUM_3DAI = members("fusarium_3dai")
BACKGROUND_3DAI = members("background_3dai")
MAPPING_3DAI = members("mapping_3dai")
CONTEXT_13DAI = members("context_13dai")
MACROPHOMINA_3DAI = members("macrophomina_3dai")
MACROPHOMINA_MAPPING = members("macrophomina_mapping_3dai")
ALL_FT_SAMPLES = sorted(set(MAPPING_3DAI + CONTEXT_13DAI))
ALL_SAMPLES = SAMPLES["sample_id"].tolist()

SAMPLE_PATTERN = "|".join(re.escape(value) for value in ALL_SAMPLES)
REFERENCE_PATTERN = "|".join(re.escape(value) for value in REFERENCES["label"])

RAW = config["outputs"]["raw"]
REF_ROOT = config["outputs"]["references"]
WORK = config["outputs"]["work"]
RESULTS = config["outputs"]["root"]
SCRIPTS = "workflow/scripts"

HOST_LABEL = config["reference_labels"]["host"]
TARGET_LABEL = config["reference_labels"]["target"]
TARGET_ALT_LABEL = config["reference_labels"]["target_alternative"]
SOURCE_LABEL = config["reference_labels"]["source_annotation"]
COMPARISON_LABELS = config["reference_labels"]["comparisons"]
MACROPHOMINA_LABEL = config["reference_labels"]["macrophomina"]

EXTENSION = int(config["analysis"]["gene_extension_bp"])
MIN_MAPQ = int(config["analysis"]["minimum_mapping_quality"])
SEED = int(config["analysis"]["random_seed"])

COMPETITIVE_FASTA = f"{WORK}/reference/competitive.fna"
CONTIG_MAP = f"{WORK}/reference/contigs.tsv"
FT_BED = f"{WORK}/reference/fusarium_thapsinum.bed"
MP_BED = f"{WORK}/reference/macrophomina_phaseolina.bed"
INDEX_PREFIX = f"{WORK}/reference/hisat2/competitive"
INDEX_FILES = [f"{INDEX_PREFIX}.{index}.ht2" for index in range(1, 9)]

FT_ALIGN = f"{WORK}/alignments/fusarium_thapsinum"
MP_ALIGN = f"{WORK}/alignments/macrophomina_phaseolina"
FT_MODELS = f"{WORK}/gene_models/fusarium_thapsinum"
MP_MODELS = f"{WORK}/gene_models/macrophomina_phaseolina"
FT_COUNTS = f"{WORK}/counts/fusarium_3dai"
FT_COUNT_MODEL = f"{RESULTS}/fusarium/count_model"
FT_GENE_QC = f"{RESULTS}/fusarium/gene_qc"
EXPRESSION_INPUTS = f"{WORK}/expression/inputs"
EXPRESSION_MODELS = f"{RESULTS}/fusarium/expression/models"
EXPRESSION_SUMMARY = f"{RESULTS}/fusarium/expression/summary"
MODULE_AUDIT = f"{RESULTS}/fusarium/module_audit"
FOCUS_ANNOTATION = f"{RESULTS}/fusarium/focus_annotation"
REGULATORY = f"{RESULTS}/fusarium/regulatory_architecture"
MOTIF_CONSERVATION = f"{RESULTS}/fusarium/motif_conservation"
CONTEXT_INPUTS = f"{WORK}/context_13dai/inputs"
CONTEXT_MODELS = f"{WORK}/context_13dai/models"
CONTEXT_RESULTS = f"{RESULTS}/fusarium/context_13dai"
HOST_COUPLING = f"{RESULTS}/host_module_coupling"
MACROPHOMINA_RESULTS = f"{RESULTS}/macrophomina/feasibility"


rule all:
    input:
        f"{RESULTS}/validation/manifests.json",
        MODULE_AUDIT,
        FOCUS_ANNOTATION,
        REGULATORY,
        MOTIF_CONSERVATION,
        CONTEXT_RESULTS,
        HOST_COUPLING,


rule macrophomina:
    input:
        MACROPHOMINA_RESULTS,


rule validate_manifests:
    input:
        samples=config["samples"],
        sets=config["analysis_sets"],
        references=config["references"],
    output:
        f"{RESULTS}/validation/manifests.json"
    shell:
        """
        python {SCRIPTS}/validate_inputs.py \
          --samples {input.samples:q} --analysis-sets {input.sets:q} \
          --references {input.references:q} --output {output:q}
        """


rule manifest_mapping_3dai:
    input:
        samples=config["samples"],
        sets=config["analysis_sets"],
    output:
        f"{WORK}/manifests/mapping_3dai.tsv"
    shell:
        """
        python {SCRIPTS}/select_samples.py \
          --samples {input.samples:q} --analysis-sets {input.sets:q} \
          --set fusarium_3dai:FUS --set background_3dai:PDB_PRIMARY \
          --output {output:q}
        """


rule manifest_context_13dai:
    input:
        samples=config["samples"],
        sets=config["analysis_sets"],
    output:
        f"{WORK}/manifests/context_13dai.tsv"
    shell:
        """
        python {SCRIPTS}/select_samples.py \
          --samples {input.samples:q} --analysis-sets {input.sets:q} \
          --set context_13dai:CONTEXT --output {output:q}
        """


rule manifest_macrophomina_3dai:
    input:
        samples=config["samples"],
        sets=config["analysis_sets"],
    output:
        f"{WORK}/manifests/macrophomina_3dai.tsv"
    shell:
        """
        python {SCRIPTS}/select_samples.py \
          --samples {input.samples:q} --analysis-sets {input.sets:q} \
          --set macrophomina_3dai:MACRO --set background_3dai:PDB_PRIMARY \
          --output {output:q}
        """


rule fetch_reference:
    output:
        genome=f"{REF_ROOT}/{{label}}/genome.fna",
        annotation=f"{REF_ROOT}/{{label}}/annotation.gff3",
        protein=f"{REF_ROOT}/{{label}}/protein.faa",
        metadata=f"{REF_ROOT}/{{label}}/metadata.json",
    params:
        accession=lambda wildcards: REFERENCE_INFO[wildcards.label]["accession"],
        require_annotation=lambda wildcards: REFERENCE_INFO[wildcards.label]["require_annotation"],
        require_protein=lambda wildcards: REFERENCE_INFO[wildcards.label]["require_protein"],
    wildcard_constraints:
        label=REFERENCE_PATTERN
    shell:
        """
        python {SCRIPTS}/fetch_reference.py \
          --accession {params.accession:q} --label {wildcards.label:q} \
          --require-annotation {params.require_annotation} \
          --require-protein {params.require_protein} \
          --output-dir {REF_ROOT}/{wildcards.label}
        """


rule competitive_reference:
    input:
        references=config["references"],
        genomes=expand(f"{REF_ROOT}/{{label}}/genome.fna", label=REFERENCES.loc[pd.to_numeric(REFERENCES["include_in_competitive"]).eq(1), "label"].tolist()),
    output:
        fasta=COMPETITIVE_FASTA,
        contigs=CONTIG_MAP,
        ft_bed=FT_BED,
        mp_bed=MP_BED,
        fai=COMPETITIVE_FASTA + ".fai",
    shell:
        """
        python {SCRIPTS}/build_competitive_reference.py \
          --references {input.references:q} --reference-root {REF_ROOT} \
          --output-fasta {output.fasta:q} --contig-map {output.contigs:q} \
          --ft-bed {output.ft_bed:q} --mp-bed {output.mp_bed:q}
        samtools faidx {output.fasta:q}
        """


rule index_competitive_reference:
    input:
        COMPETITIVE_FASTA
    output:
        INDEX_FILES
    threads:
        config["threads"]["reference_index"]
    params:
        prefix=INDEX_PREFIX
    shell:
        "hisat2-build -p {threads} {input:q} {params.prefix:q}"


rule download_fastq:
    input:
        metadata=config["samples"]
    output:
        fastq=f"{RAW}/{{sample}}.fastq.gz",
        metrics=f"{RESULTS}/downloads/{{sample}}.tsv",
    params:
        url=lambda wildcards: SAMPLE_INFO[wildcards.sample]["fastq_url"],
        md5=lambda wildcards: SAMPLE_INFO[wildcards.sample]["fastq_md5"],
        reads=lambda wildcards: SAMPLE_INFO[wildcards.sample]["expected_reads"],
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    shell:
        """
        python {SCRIPTS}/download_fastq.py \
          --sample {wildcards.sample:q} --url {params.url:q} --md5 {params.md5:q} \
          --expected-reads {params.reads} --output {output.fastq:q} --metrics {output.metrics:q}
        """


rule align_fusarium_thapsinum:
    input:
        fastq=f"{RAW}/{{sample}}.fastq.gz",
        index=INDEX_FILES,
        contigs=CONTIG_MAP,
        bed=FT_BED,
    output:
        bam=f"{FT_ALIGN}/bam/{{sample}}.bam",
        bai=f"{FT_ALIGN}/bam/{{sample}}.bam.bai",
        species=f"{FT_ALIGN}/species/{{sample}}.tsv",
        summary=f"{FT_ALIGN}/summary/{{sample}}.txt",
        fastp=f"{FT_ALIGN}/fastp/{{sample}}.json",
        coverage=f"{FT_ALIGN}/coverage/{{sample}}.thresholds.bed.gz",
    threads:
        config["threads"]["align"]
    wildcard_constraints:
        sample="|".join(re.escape(value) for value in ALL_FT_SAMPLES)
    shell:
        """
        bash {SCRIPTS}/align_sample.sh \
          {input.fastq:q} {wildcards.sample:q} {INDEX_PREFIX:q} {input.contigs:q} {input.bed:q} \
          FT {MIN_MAPQ} {threads} {output.bam:q} {output.species:q} {output.summary:q} \
          {output.fastp:q} {output.coverage:q}
        """


rule align_macrophomina_phaseolina:
    input:
        fastq=f"{RAW}/{{sample}}.fastq.gz",
        index=INDEX_FILES,
        contigs=CONTIG_MAP,
        bed=MP_BED,
    output:
        bam=f"{MP_ALIGN}/bam/{{sample}}.bam",
        bai=f"{MP_ALIGN}/bam/{{sample}}.bam.bai",
        species=f"{MP_ALIGN}/species/{{sample}}.tsv",
        summary=f"{MP_ALIGN}/summary/{{sample}}.txt",
        fastp=f"{MP_ALIGN}/fastp/{{sample}}.json",
        coverage=f"{MP_ALIGN}/coverage/{{sample}}.thresholds.bed.gz",
    threads:
        config["threads"]["align"]
    wildcard_constraints:
        sample="|".join(re.escape(value) for value in MACROPHOMINA_MAPPING)
    shell:
        """
        bash {SCRIPTS}/align_sample.sh \
          {input.fastq:q} {wildcards.sample:q} {INDEX_PREFIX:q} {input.contigs:q} {input.bed:q} \
          MP {MIN_MAPQ} {threads} {output.bam:q} {output.species:q} {output.summary:q} \
          {output.fastp:q} {output.coverage:q}
        """


rule summarize_mapping_3dai:
    input:
        manifest=rules.manifest_mapping_3dai.output,
        species=expand(f"{FT_ALIGN}/species/{{sample}}.tsv", sample=MAPPING_3DAI),
        summaries=expand(f"{FT_ALIGN}/summary/{{sample}}.txt", sample=MAPPING_3DAI),
        fastp=expand(f"{FT_ALIGN}/fastp/{{sample}}.json", sample=MAPPING_3DAI),
        coverage=expand(f"{FT_ALIGN}/coverage/{{sample}}.thresholds.bed.gz", sample=MAPPING_3DAI),
    output:
        directory(f"{RESULTS}/fusarium/mapping_3dai")
    shell:
        """
        python {SCRIPTS}/summarize_mapping.py \
          --manifest {input.manifest:q} --species-dir {FT_ALIGN}/species \
          --summary-dir {FT_ALIGN}/summary --fastp-dir {FT_ALIGN}/fastp \
          --coverage-dir {FT_ALIGN}/coverage --target-group FT --inoculated-label fus \
          --rpm-basis input --output-dir {output:q}
        """


rule summarize_context_mapping:
    input:
        manifest=rules.manifest_context_13dai.output,
        species=expand(f"{FT_ALIGN}/species/{{sample}}.tsv", sample=CONTEXT_13DAI),
        summaries=expand(f"{FT_ALIGN}/summary/{{sample}}.txt", sample=CONTEXT_13DAI),
        fastp=expand(f"{FT_ALIGN}/fastp/{{sample}}.json", sample=CONTEXT_13DAI),
        coverage=expand(f"{FT_ALIGN}/coverage/{{sample}}.thresholds.bed.gz", sample=CONTEXT_13DAI),
    output:
        directory(f"{WORK}/context_13dai/mapping")
    shell:
        """
        python {SCRIPTS}/summarize_mapping.py \
          --manifest {input.manifest:q} --species-dir {FT_ALIGN}/species \
          --summary-dir {FT_ALIGN}/summary --fastp-dir {FT_ALIGN}/fastp \
          --coverage-dir {FT_ALIGN}/coverage --target-group FT --inoculated-label fus \
          --rpm-basis post_filter --output-dir {output:q}
        """


rule summarize_macrophomina_mapping:
    input:
        manifest=rules.manifest_macrophomina_3dai.output,
        species=expand(f"{MP_ALIGN}/species/{{sample}}.tsv", sample=MACROPHOMINA_MAPPING),
        summaries=expand(f"{MP_ALIGN}/summary/{{sample}}.txt", sample=MACROPHOMINA_MAPPING),
        fastp=expand(f"{MP_ALIGN}/fastp/{{sample}}.json", sample=MACROPHOMINA_MAPPING),
        coverage=expand(f"{MP_ALIGN}/coverage/{{sample}}.thresholds.bed.gz", sample=MACROPHOMINA_MAPPING),
    output:
        directory(f"{WORK}/macrophomina/mapping")
    shell:
        """
        python {SCRIPTS}/summarize_mapping.py \
          --manifest {input.manifest:q} --species-dir {MP_ALIGN}/species \
          --summary-dir {MP_ALIGN}/summary --fastp-dir {MP_ALIGN}/fastp \
          --coverage-dir {MP_ALIGN}/coverage --target-group MP --inoculated-label macro \
          --rpm-basis post_filter --output-dir {output:q}
        """


rule fusarium_gene_models:
    input:
        genome=f"{REF_ROOT}/{TARGET_LABEL}/genome.fna",
        source_proteins=f"{REF_ROOT}/{SOURCE_LABEL}/protein.faa",
        fai=COMPETITIVE_FASTA + ".fai",
    output:
        directory(FT_MODELS)
    threads:
        config["threads"]["annotation"]
    shell:
        """
        bash {SCRIPTS}/project_gene_models.sh \
          {input.genome:q} {input.source_proteins:q} {input.fai:q} \
          {TARGET_LABEL}__ target_genes {EXTENSION} {threads} {output:q}
        """


rule macrophomina_gene_models:
    input:
        annotation=f"{REF_ROOT}/{MACROPHOMINA_LABEL}/annotation.gff3",
        fai=COMPETITIVE_FASTA + ".fai",
    output:
        directory(MP_MODELS)
    shell:
        """
        python {SCRIPTS}/build_quantseq_saf.py \
          --gff {input.annotation:q} --fai {input.fai:q} \
          --prefix {MACROPHOMINA_LABEL}__ --already-prefixed 0 --extension {EXTENSION} \
          --output-prefix target_genes --output-dir {output:q}
        """


rule count_fusarium_genes:
    input:
        manifest=rules.manifest_mapping_3dai.output,
        bams=expand(f"{FT_ALIGN}/bam/{{sample}}.bam", sample=MAPPING_3DAI),
        models=rules.fusarium_gene_models.output,
    output:
        directory(FT_COUNTS)
    threads:
        config["threads"]["feature_counts"]
    shell:
        """
        bash {SCRIPTS}/count_genes.sh \
          {input.manifest:q} {FT_ALIGN}/bam {input.models:q} target_genes \
          {EXTENSION} {threads} {output:q}
        """


rule select_fusarium_count_model:
    input:
        counts=rules.count_fusarium_genes.output,
        manifest=rules.manifest_mapping_3dai.output,
    output:
        directory(FT_COUNT_MODEL)
    shell:
        """
        python {SCRIPTS}/select_count_model.py \
          --counts-dir {input.counts:q} --manifest {input.manifest:q} \
          --extension {EXTENSION} --inoculated-label fus --control-set-label PDB_PRIMARY \
          --output-dir {output:q}
        """


rule summarize_fusarium_gene_counts:
    input:
        manifest=rules.manifest_mapping_3dai.output,
        selection=rules.select_fusarium_count_model.output,
        metadata=rules.fusarium_gene_models.output,
    output:
        directory(FT_GENE_QC)
    shell:
        """
        python {SCRIPTS}/summarize_gene_counts.py \
          --manifest {input.manifest:q} \
          --counts {input.selection}/selected_counts.tsv \
          --sample-metrics {input.selection}/selected_sample_metrics.tsv \
          --gene-metadata {input.metadata}/target_genes_metadata.tsv \
          --min-count {config[analysis][minimum_gene_count]} \
          --min-samples {config[analysis][minimum_gene_samples]} \
          --inoculated-label fus --control-label pdb --output-dir {output:q}
        """


rule prepare_expression_inputs:
    input:
        selection=rules.select_fusarium_count_model.output,
        manifest=rules.manifest_mapping_3dai.output,
        qc=rules.summarize_fusarium_gene_counts.output,
        metadata=rules.fusarium_gene_models.output,
        mapping=rules.summarize_mapping_3dai.output,
        source_proteins=f"{REF_ROOT}/{SOURCE_LABEL}/protein.faa",
        terms=config["functional_terms"],
    output:
        directory(EXPRESSION_INPUTS)
    shell:
        """
        python {SCRIPTS}/prepare_expression_inputs.py \
          --counts {input.selection}/selected_counts.tsv --manifest {input.manifest:q} \
          --sample-qc {input.qc}/sample_gene_qc.tsv --mapping {input.mapping}/sample_mapping.tsv \
          --background {input.qc}/background_genes.tsv --usable-genes {input.qc}/usable_genes.tsv \
          --gene-metadata {input.metadata}/target_genes_metadata.tsv \
          --source-proteins {input.source_proteins:q} --functional-terms {input.terms:q} \
          --min-count {config[analysis][minimum_gene_count]} \
          --min-samples {config[analysis][minimum_gene_samples]} --output-dir {output:q}
        """


rule fit_expression_models:
    input:
        data=rules.prepare_expression_inputs.output,
        mapping=rules.summarize_mapping_3dai.output,
    output:
        directory(EXPRESSION_MODELS)
    shell:
        """
        Rscript {SCRIPTS}/fit_expression_models.R \
          {input.data}/counts_primary.tsv {input.data}/samples_primary.tsv \
          {input.data}/counts_all.tsv {input.data}/samples_all.tsv \
          {input.data}/functional_gene_sets.tsv {input.mapping}/sample_mapping.tsv \
          {output:q} {config[analysis][false_discovery_rate]} \
          {config[analysis][minimum_absolute_log2_fold_change]}
        """


rule summarize_expression_models:
    input:
        models=rules.fit_expression_models.output,
        data=rules.prepare_expression_inputs.output,
    output:
        directory(EXPRESSION_SUMMARY)
    shell:
        """
        python {SCRIPTS}/summarize_expression.py \
          --model-dir {input.models:q} --annotation {input.data}/gene_annotation.tsv \
          --gene-sets {input.data}/functional_gene_sets.tsv --output-dir {output:q} \
          --fdr {config[analysis][false_discovery_rate]} \
          --min-abs-log2fc {config[analysis][minimum_absolute_log2_fold_change]}
        """


rule audit_aromatic_module:
    input:
        data=rules.prepare_expression_inputs.output,
        models=rules.fit_expression_models.output,
        summary=rules.summarize_expression_models.output,
    output:
        directory(MODULE_AUDIT)
    shell:
        """
        python {SCRIPTS}/audit_aromatic_module.py \
          --samples-primary {input.data}/samples_primary.tsv \
          --samples-all {input.data}/samples_all.tsv \
          --counts {input.data}/counts_primary.tsv \
          --logcpm-primary {input.models}/primary_state/logCPM.tsv \
          --logcpm-all {input.models}/all_state/logCPM.tsv \
          --focus-genes {input.summary}/focus_genes.tsv \
          --de-contrast {input.summary}/contrasts/bmr6_vs_wt_wet.tsv \
          --focus-module aromatic_phenolic_metabolism \
          --seed {SEED} --matched-iterations {config[analysis][matched_set_iterations]} \
          --output-dir {output:q}
        """


rule annotate_focus_proteins:
    input:
        module=rules.audit_aromatic_module.output,
        data=rules.prepare_expression_inputs.output,
        genome=f"{REF_ROOT}/{TARGET_LABEL}/genome.fna",
        comparison_proteins=expand(f"{REF_ROOT}/{{label}}/protein.faa", label=[label for label in COMPARISON_LABELS if label != TARGET_ALT_LABEL]),
    output:
        directory(FOCUS_ANNOTATION)
    threads:
        config["threads"]["annotation"]
    params:
        comparisons=" ".join(
            f"--comparison {label}={REF_ROOT}/{label}/protein.faa"
            for label in COMPARISON_LABELS if label != TARGET_ALT_LABEL
        )
    shell:
        """
        python {SCRIPTS}/annotate_focus_proteins.py \
          --focus-genes {EXPRESSION_SUMMARY}/focus_genes.tsv \
          --source-proteins {input.data}/projected_source_proteins.faa \
          --fl4-genome {input.genome:q} {params.comparisons} \
          --work-dir {WORK}/focus_annotation --threads {threads} --output-dir {output:q}
        """


rule regulatory_architecture:
    input:
        code=f"{SCRIPTS}/analyze_regulatory_architecture.py",
        motif_code=f"{SCRIPTS}/motif_inference.py",
        data=rules.prepare_expression_inputs.output,
        models=rules.fit_expression_models.output,
        summary=rules.summarize_expression_models.output,
        module=rules.audit_aromatic_module.output,
        genome=f"{REF_ROOT}/{TARGET_LABEL}/genome.fna",
    output:
        directory(REGULATORY)
    params:
        motif_lengths=" ".join(str(k) for k in config["analysis"]["motif_lengths"]),
        module_name=config["analysis"]["focus_module"]
    shell:
        """
        python {SCRIPTS}/analyze_regulatory_architecture.py \
          --annotation {input.data}/gene_annotation.tsv \
          --focus-genes {input.summary}/focus_genes.tsv \
          --module-scores {input.module}/module_scores.tsv \
          --logcpm {input.models}/primary_state/logCPM.tsv \
          --counts {input.data}/counts_primary.tsv \
          --samples {input.data}/samples_primary.tsv \
          --interaction-de {input.summary}/contrasts/bmr6_water_interaction.tsv \
          --genome {input.genome:q} --seed {SEED} \
          --module-name {params.module_name:q} --motif-lengths {params.motif_lengths} \
          --motif-fdr {config[analysis][motif_false_discovery_rate]} \
          --matched-iterations {config[analysis][regulatory_matched_set_iterations]} \
          --regulator-permutations {config[analysis][regulator_permutations]} \
          --promoter-bp {config[analysis][promoter_length_bp]} --output-dir {output:q}
        """


rule motif_conservation:
    input:
        code=f"{SCRIPTS}/test_motif_conservation.py",
        motif_code=f"{SCRIPTS}/motif_inference.py",
        regulatory=rules.regulatory_architecture.output,
        data=rules.prepare_expression_inputs.output,
        genomes=expand(f"{REF_ROOT}/{{label}}/genome.fna", label=COMPARISON_LABELS),
    output:
        directory(MOTIF_CONSERVATION)
    threads:
        config["threads"]["annotation"]
    params:
        genomes=" ".join(f"--genome {label}={REF_ROOT}/{label}/genome.fna" for label in COMPARISON_LABELS),
        mode=config["analysis"]["motif_follow_up"]["mode"],
        motif=config["analysis"]["motif_follow_up"]["motif"]
    shell:
        """
        python {SCRIPTS}/test_motif_conservation.py \
          --motif-table {input.regulatory}/tables/promoter_kmer_enrichment.tsv \
          --selection-mode {params.mode:q} --motif {params.motif:q} \
          --selection-fdr {config[analysis][motif_false_discovery_rate]} \
          --reference-label {TARGET_LABEL:q} --within-species-label {TARGET_ALT_LABEL:q} \
          --focus-promoters {input.regulatory}/tables/focus_promoters.tsv \
          --matched-pools {input.regulatory}/tables/matched_background_pools.tsv \
          --focus-promoter-fasta {input.regulatory}/sequences/focus_promoters.fasta \
          --background-promoter-fasta {input.regulatory}/sequences/matched_background_promoters.fasta \
          --source-proteins {input.data}/projected_source_proteins.faa \
          --gene-annotation {input.data}/gene_annotation.tsv {params.genomes} \
          --work-dir {WORK}/motif_conservation --threads {threads} \
          --promoter-bp {config[analysis][promoter_length_bp]} \
          --matched-iterations {config[analysis][motif_conservation_iterations]} \
          --seed {SEED} --output-dir {output:q}
        """


rule count_context_genes:
    input:
        manifest=rules.manifest_context_13dai.output,
        bams=expand(f"{FT_ALIGN}/bam/{{sample}}.bam", sample=CONTEXT_13DAI),
        models=rules.fusarium_gene_models.output,
        selection=rules.select_fusarium_count_model.output,
    output:
        counts=f"{WORK}/context_13dai/counts.txt",
        summary=f"{WORK}/context_13dai/counts.txt.summary",
    threads:
        config["threads"]["feature_counts"]
    shell:
        """
        python {SCRIPTS}/count_selected_genes.py \
          --manifest {input.manifest:q} --bam-dir {FT_ALIGN}/bam \
          --gene-model-dir {input.models:q} --model-prefix target_genes \
          --selection {input.selection}/count_selection.json \
          --threads {threads} --output {output.counts:q}
        """


rule prepare_context_13dai:
    input:
        manifest=rules.manifest_context_13dai.output,
        mapping=rules.summarize_context_mapping.output,
        counts=rules.count_context_genes.output.counts,
        summary=rules.count_context_genes.output.summary,
        focus=rules.summarize_expression_models.output,
        prior=rules.summarize_fusarium_gene_counts.output,
    output:
        directory(CONTEXT_INPUTS)
    shell:
        """
        python {SCRIPTS}/prepare_13dai_context.py \
          --manifest {input.manifest:q} --mapping {input.mapping}/sample_mapping.tsv \
          --counts {input.counts:q} --featurecounts-summary {input.summary:q} \
          --focus-genes {input.focus}/focus_genes.tsv \
          --prior-background {input.prior}/background_genes.tsv \
          --minimum-assigned-reads {config[analysis][context_minimum_assigned_reads]} \
          --minimum-detected-genes {config[analysis][context_minimum_detected_genes]} \
          --minimum-target-reads {config[analysis][context_minimum_target_reads]} \
          --minimum-specificity {config[analysis][context_minimum_target_specificity]} \
          --output-dir {output:q}
        """


rule fit_context_models:
    input:
        data=rules.prepare_context_13dai.output,
    output:
        directory(CONTEXT_MODELS)
    shell:
        """
        Rscript {SCRIPTS}/fit_context_models.R \
          {input.data}/counts_primary.tsv {input.data}/samples_primary.tsv \
          {input.data}/gene_filter.tsv {input.data}/focus_genes.tsv \
          {output:q} {config[analysis][minimum_gene_count]}
        """


rule analyze_context_13dai:
    input:
        data=rules.prepare_context_13dai.output,
        models=rules.fit_context_models.output,
    output:
        directory(CONTEXT_RESULTS)
    shell:
        """
        python {SCRIPTS}/analyze_13dai_context.py \
          --sample-summary {input.data}/sample_summary.tsv \
          --state-logcpm {input.models}/state/logCPM.tsv \
          --focus-genes {input.data}/focus_genes.tsv \
          --state-camera {input.models}/state/camera.tsv \
          --tissue-camera {input.models}/tissue/camera.tsv \
          --output-dir {output:q}
        """


rule fetch_host_eigengenes:
    output:
        f"{WORK}/host_eigengenes.csv"
    params:
        url=config["host_eigengene_url"]
    shell:
        """
        python {SCRIPTS}/fetch_file.py --url {params.url:q} --minimum-bytes 1000 --output {output:q}
        """


rule host_module_coupling:
    input:
        eigengenes=rules.fetch_host_eigengenes.output,
        module=rules.audit_aromatic_module.output,
        models=rules.fit_expression_models.output,
        data=rules.prepare_expression_inputs.output,
        traits=config["samples"],
        focus=rules.summarize_expression_models.output,
    output:
        directory(HOST_COUPLING)
    shell:
        """
        python {SCRIPTS}/host_module_coupling.py \
          --eigengenes {input.eigengenes:q} \
          --module-scores {input.module}/module_scores.tsv \
          --focus-genes {input.focus}/focus_genes.tsv \
          --logcpm-all {input.models}/all_state/logCPM.tsv \
          --samples-primary {input.data}/samples_primary.tsv \
          --samples-all {input.data}/samples_all.tsv \
          --sample-traits {input.traits:q} \
          --permutations {config[analysis][host_module_permutations]} \
          --seed {SEED} --output-dir {output:q}
        """


rule count_macrophomina_genes:
    input:
        manifest=rules.manifest_macrophomina_3dai.output,
        bams=expand(f"{MP_ALIGN}/bam/{{sample}}.bam", sample=MACROPHOMINA_MAPPING),
        models=rules.macrophomina_gene_models.output,
    output:
        directory(f"{WORK}/counts/macrophomina_3dai")
    threads:
        config["threads"]["feature_counts"]
    shell:
        """
        bash {SCRIPTS}/count_genes.sh \
          {input.manifest:q} {MP_ALIGN}/bam {input.models:q} target_genes \
          {EXTENSION} {threads} {output:q}
        """


rule select_macrophomina_count_model:
    input:
        counts=rules.count_macrophomina_genes.output,
        manifest=rules.manifest_macrophomina_3dai.output,
    output:
        directory(f"{WORK}/macrophomina/count_model")
    shell:
        """
        python {SCRIPTS}/select_count_model.py \
          --counts-dir {input.counts:q} --manifest {input.manifest:q} \
          --extension {EXTENSION} --inoculated-label macro --control-set-label PDB_PRIMARY \
          --output-dir {output:q}
        """


rule assess_macrophomina_signal:
    input:
        manifest=rules.manifest_macrophomina_3dai.output,
        mapping=rules.summarize_macrophomina_mapping.output,
        selection=rules.select_macrophomina_count_model.output,
    output:
        directory(MACROPHOMINA_RESULTS)
    shell:
        """
        python {SCRIPTS}/assess_target_signal.py \
          --manifest {input.manifest:q} --mapping {input.mapping}/sample_mapping.tsv \
          --counts {input.selection}/selected_counts.tsv \
          --featurecounts-summary {input.selection}/selected_featurecounts_summary.tsv \
          --inoculated-label macro --control-label pdb \
          --bam-suffix .bam --output-dir {output:q}
        """


include: "workflow/rules/sample_sensitivity.smk"
