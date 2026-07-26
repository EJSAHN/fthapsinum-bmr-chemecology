# fthapsinum-bmr-chemecology

Reproducible pathogen-side RNA-seq workflow accompanying the manuscript:

**Water-dependent deployment of a dispersed aromatic-processing programme in *Fusarium thapsinum* during brown-midrib sorghum colonisation**

The workflow recovers and analyses *F. thapsinum* RNA from a public sorghum QuantSeq experiment and tests whether host genotype and water treatment are associated with fungal RNA burden, relative transcript allocation, aromatic- and phenolic-processing genes, promoter architecture and conservation of a candidate promoter element.

The repository remains private during USDA institutional clearance and will be released publicly before publication.

## Public data

The experiment includes wild-type Tx430, `bmr6` and `bmr12` sorghum under water-limited and well-watered conditions. Sequence data are available under NCBI BioProject **PRJNA573931**. The source publication is identified by DOI **10.1186/s12870-021-03149-5**.

## Workflow

The default workflow:

1. validates the sample and reference manifests;
2. downloads public FASTQ files and genome assemblies;
3. builds a prefixed competitive reference containing sorghum, *F. thapsinum*, *Macrophomina phaseolina* and close *Fusarium* references;
4. aligns single-end QuantSeq reads and summarises species-resolved primary alignments;
5. projects *F. verticillioides* proteins onto the *F. thapsinum* FL-4 genome;
6. compares native and 3′-extended gene models across featureCounts strand settings;
7. filters recurrent control-associated loci and low-information libraries;
8. fits robust edgeR quasi-likelihood models on fungal-state and tissue-output scales;
9. audits the aromatic/phenolic-processing gene set using exact permutation, burden-adjusted, leave-one-out and matched-set tests;
10. evaluates protein annotation, genomic dispersion, promoter enrichment, candidate regulators and motif conservation;
11. tests the fungal gene set in the 13-day context; and
12. tests associations with published host WGCNA eigengenes.

An optional target evaluates whether the *M. phaseolina* libraries contain enough reference-supported pathogen RNA for gene-level analysis. It reports feasibility and does not force a differential-expression model when pathogen counts are insufficient.

## Installation

Create the software environment with Mamba or Micromamba:

```bash
mamba env create -f environment.yml
mamba activate fthapsinum-bmr-chemecology
```

Run the complete *F. thapsinum* analysis from the repository root:

```bash
snakemake --cores 12
```

Run the optional *M. phaseolina* feasibility analysis:

```bash
snakemake macrophomina --cores 12
```

A dry run shows the planned jobs without downloading or computing data:

```bash
snakemake --dry-run --printshellcmds
```

All paths are relative to the repository root. Thread counts, thresholds, accessions, URLs and analysis parameters are stored in `config.yaml` and the files under `resources/`.

## Repository layout

- `Snakefile`: workflow graph and analysis targets
- `config.yaml`: analysis parameters and reference roles
- `resources/`: sample, analysis-set, reference and functional-term manifests
- `workflow/scripts/`: command-line analysis programs
- `tests/`: manifest, helper, portability and policy tests
- `.github/workflows/ci.yml`: syntax checks, Snakemake dry run and tests

## Inputs

### Public sequence data

`resources/samples.tsv` records the public sample name, canonical sample name, experimental factors, SRA run, ENA FASTQ URL, MD5 checksum, expected read count and available traits for 129 libraries.

Three published `bmr6` PDB/well-watered sample names omit the day token. Their canonical identifiers and 3-day assignment are recorded explicitly with `timepoint_basis=inferred_from_published_notebook`. The primary pathogen-background set uses the 21 controls whose day is stated directly in the sample identifier. The inferred labels are retained for the balanced host-eigengene analysis.

Twenty-eight libraries, including 0-day baselines and additional PDB replicates, are retained in the complete public manifest for provenance but are not assigned to a downstream analysis set.

### Analysis sets

`resources/analysis_sets.tsv` defines all sample selections used by the workflow. Selection logic is kept out of the analysis scripts.

| Set | Libraries | Purpose |
| --- | ---: | --- |
| `fusarium_3dai` | 24 | six genotype–water cells with four inoculated libraries each |
| `background_3dai` | 21 | primary PDB background panel |
| `mapping_3dai` | 45 | competitive mapping and fungal gene-count selection |
| `context_13dai` | 29 | 15 inoculated and 14 PDB libraries at 13 days |
| `host_3dai` | 48 | balanced FUS/PDB host-eigengene design |
| `macrophomina_3dai` | 24 | optional pathogen-side feasibility analysis |
| `macrophomina_mapping_3dai` | 45 | 24 inoculated libraries plus the 21 primary PDB controls |

### References

Genome accessions and their roles are listed in `resources/references.tsv`. The workflow downloads assemblies with NCBI Datasets and records the selected genomic FASTA, annotation, protein FASTA and genome size.

## Outputs

The workflow writes tables, JSON summaries, FASTA files, model objects and software-session records. It does not generate manuscript figures.

Primary output locations are:

```text
results/fusarium/mapping_3dai/
results/fusarium/count_model/
results/fusarium/gene_qc/
results/fusarium/expression/
results/fusarium/module_audit/
results/fusarium/focus_annotation/
results/fusarium/regulatory_architecture/
results/fusarium/motif_conservation/
results/fusarium/context_13dai/
results/host_module_coupling/
```

The optional *M. phaseolina* feasibility result is written to:

```text
results/macrophomina/feasibility/
```

Large FASTQ, reference, index and BAM files are excluded by `.gitignore`.

## Statistical design

The primary fungal expression analysis uses 22 quality-passing inoculated libraries. The 24-library set is retained as a sensitivity analysis. Fungal-state models use TMMwsp normalisation within fungal counts. Tissue-output models use total post-filtering tissue reads as the library-size offset. Contrasts are estimated with robust edgeR quasi-likelihood models.

The aromatic/phenolic-processing gene set is evaluated as a fixed set after its definition. Exact tests preserve the factorial exchangeability structure. Additional checks include fungal-burden adjustment, leave-one-sample-out analysis, leave-one-gene-out analysis and expression- and detection-matched random gene sets.

Promoter analysis uses gene-oriented 1-kb upstream sequences. K-mer enrichment is evaluated against matched background promoters. Motif conservation is tested in a second *F. thapsinum* assembly and three related *Fusarium* genomes. Enrichment identifies a candidate cis-element; it does not establish transcription-factor binding.

## Interpretation limits

- RNA-derived pathogen signal is not a direct biomass measurement.
- The main expression result is gene-set level; it is not supported by individually significant fungal-state genes at the stated threshold.
- Association with `bmr6` under well-watered conditions does not identify the inducing host compound.
- The promoter k-mer is a candidate regulatory element, not a validated binding site.
- The public *M. phaseolina* libraries contain too little reference-supported pathogen RNA for a cross-pathogen expression comparison.
- Published host eigengenes do not show a robust adjusted association with the fungal aromatic/phenolic-processing score.

## Tests

Run repository tests with:

```bash
pytest -q
```

The test suite checks manifests, analysis-set membership, helper functions, script syntax and repository policy constraints. The GitHub Actions workflow also builds the Snakemake dry-run DAG.

## Citation

When using this workflow, cite the accompanying manuscript, NCBI BioProject **PRJNA573931**, the source experiment (**10.1186/s12870-021-03149-5**), the genome accessions in `resources/references.tsv` and the software packages reported in the workflow outputs.

A machine-readable software citation is provided in `CITATION.cff`.

## Licence

Code is released under the MIT License. Public sequence data and genome assemblies remain subject to the terms of their source repositories.
