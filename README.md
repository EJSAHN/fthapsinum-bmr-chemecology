# fthapsinum-bmr-chemecology

Pathogen-side RNA-seq analysis of *Fusarium thapsinum* during sorghum stalk colonisation.

Associated study: **Water-dependent deployment of an aromatic-processing programme in *Fusarium thapsinum* during brown-midrib sorghum colonisation**.

The workflow analyses a public sorghum QuantSeq experiment to distinguish RNA-derived fungal burden from relative fungal transcript allocation. It tests host-genotype and water-treatment contrasts, evaluates an aromatic/phenolic-processing gene set, and examines projected gene locations and promoter sequences. Code is publicly available; raw sequencing data remain in their source archive.

## Public data

NCBI BioProject **PRJNA573931** contains the source experiment described by Khasin et al. (2021), DOI **10.1186/s12870-021-03149-5**. The experiment includes Tx430 wild type, CAD-deficient `bmr6` and COMT-deficient `bmr12` under water-limited and well-watered conditions.

`resources/samples.tsv` records 129 libraries with public and canonical sample names, experimental factors, SRA accessions, ENA FASTQ URLs, MD5 checksums, expected read counts and available traits. Three published `bmr6` PDB/well-watered names omit the day token; their 3-day assignment is explicitly marked `timepoint_basis=inferred_from_published_notebook`. They enter the balanced host-eigengene analysis, not the primary 21-library background panel. Twenty-eight libraries are retained for provenance without assignment to a downstream analysis set.

`resources/analysis_sets.tsv` defines sample membership:

| Set | Libraries | Purpose |
| --- | ---: | --- |
| `fusarium_3dai` | 24 | Three genotypes, two water treatments and four infected replicates per cell |
| `background_3dai` | 21 | Run-matched PDB background panel |
| `mapping_3dai` | 45 | Competitive mapping and fungal gene-count selection |
| `context_13dai` | 29 | Fifteen infected and fourteen PDB libraries at 13 days |
| `host_3dai` | 48 | Balanced FUS/PDB host-eigengene design |
| `macrophomina_3dai` | 24 | Pathogen-side feasibility analysis |
| `macrophomina_mapping_3dai` | 45 | Twenty-four inoculated libraries plus the primary PDB controls |

Genome accessions and analysis roles are recorded in `resources/references.tsv`. Functional terms are recorded in `resources/functional_terms.yaml`.

## Installation and execution

Create and activate the environment from the repository root:

```bash
mamba env create -f environment.yml
mamba activate fthapsinum-bmr-chemecology
```

Inspect the planned jobs, then run the main workflow:

```bash
snakemake --dry-run --cores 1
snakemake --cores 12
```

The optional *Macrophomina phaseolina* feasibility target is:

```bash
snakemake macrophomina --cores 12
```

It reports insufficient pathogen RNA without forcing an expression model. Paths are relative to the repository root. Parameters, thread counts and reference roles are specified in `config.yaml` and `resources/`.

## Analyses and outputs

The workflow downloads public inputs, builds a prefixed competitive reference, aligns single-end QuantSeq reads, and summarises species-resolved primary alignments. It projects source proteins to FL-4, compares native and 3-prime-extended counting models, filters background-associated loci and low-information libraries, and fits robust edgeR quasi-likelihood models. Gene-set analyses include CAMERA, exact permutation, burden adjustment, leave-one-out checks and matched gene sets.

Further outputs cover protein annotation, sequence-record-level gene locations, residual coexpression, candidate regulators, promoter k-mers, comparative promoter observations, the 13-day context and published host eigengenes.

| Output directory | Contents |
| --- | --- |
| `results/fusarium/mapping_3dai/` | Competitive-mapping summaries |
| `results/fusarium/count_model/` | QuantSeq counting-model selection |
| `results/fusarium/gene_qc/` | Gene and library quality metrics |
| `results/fusarium/expression/` | Fungal-state and tissue-output models |
| `results/fusarium/module_audit/` | Targeted gene-set sensitivity analyses |
| `results/fusarium/focus_annotation/` | Protein-annotation evidence |
| `results/fusarium/regulatory_architecture/` | Gene locations, coexpression and complete-family motif tests |
| `results/fusarium/motif_conservation/` | Explicitly labelled motif follow-up |
| `results/fusarium/context_13dai/` | Contextual 13-day analysis |
| `results/host_module_coupling/` | Host-eigengene associations |
| `results/macrophomina/feasibility/` | Optional pathogen-signal feasibility assessment |

Outputs are tables, JSON summaries, FASTA files, model objects and software-session records. Manuscript figures are not generated. Downloaded reads, genomes, indexes, BAM files and local results are excluded by `.gitignore`.

## Statistical interpretation

The primary expression analysis uses 22 quality-passing infected libraries, with all 24 retained for sensitivity analysis. Fungal-state models use TMMwsp normalisation within fungal counts; tissue-output models use total tissue-read exposure. The ten-gene aromatic/phenolic-processing set was selected from the same experiment used for targeted follow-up. Its sensitivity analyses provide internal evidence, not independent biological replication.

### Promoter discovery and exploratory follow-up

Gene-oriented 1-kb upstream sequences are screened against expression-, detection-, GC- and length-matched background promoters. All **10,272 reverse-complement-collapsed A/C/G/T 6- and 7-mers** enter one Benjamini-Hochberg family, including words absent from the target promoters. Enrichment uses promoter presence and a one-sided Fisher test; there is no target-frequency filter before adjustment.

For the analysed target and background sequences, no motif passes the configured full-family threshold of 0.10. ATGCTAA/TTAGCAT occurs in 6/10 target promoters and 38/557 background promoters, with unadjusted Fisher P = 2.59645 x 10^-5 and full-family q = 0.266708. It is retained only as an **exploratory candidate**.

The supplied configuration explicitly requests that descriptive follow-up:

```yaml
motif_follow_up:
  mode: exploratory
  motif: ATGCTAA
```

For discovery-driven selection, set both values as follows under `analysis` in `config.yaml`:

```yaml
motif_follow_up:
  mode: auto
  motif: AUTO
```

AUTO does not substitute a failed candidate when no motif passes. It returns `NO_DISCOVERY_MOTIF` with no selected motif; multiple passing motifs also remain unselected rather than being reduced to the best-looking hit. Exploratory output carries the candidate's full-family q value and does not promote fixed-candidate empirical probabilities to search-wide significance. Comparative retention in related genomes is not an independent infection-expression experiment or a binding assay.

Details of the hypothesis family, selection modes and location summaries are in [Promoter analysis](docs/promoter_analysis.md).

### Other limits

- RNA-derived pathogen signal is not a direct biomass measurement, and adjustment for it does not exclude every abundance-related explanation.
- The primary expression result is gene-set level, without individually significant fungal-state genes at the stated threshold.
- Association with well-watered `bmr6` does not identify an inducing host compound or establish a virulence function.
- Separate scaffold records do not resolve chromosome-level distances. Location outputs distinguish measured same-record gaps from unknown cross-record distances; projected corresponding loci are not a new synteny validation.
- The *M. phaseolina* libraries do not support a cross-pathogen expression comparison. The 13-day arm lacks `bmr6`, and published host eigengenes show no robust adjusted coupling with the aromatic score.

## Tests

```bash
python -m compileall -q workflow/scripts
pytest -q
snakemake --dry-run --cores 1
snakemake macrophomina --dry-run --cores 1
```

Tests cover manifests, helpers, portability, complete-family inference, candidate-selection behaviour and no-hit output handling. GitHub Actions runs syntax checks, unit tests and both workflow dry runs. These checks do not reprocess the sequencing data or establish biological validity. The affected production programs have also been executed against the archived study inputs, including regenerated promoters and both follow-up modes.

## Citation and licence

Software authors are listed in `CITATION.cff`. Cite the source experiment, BioProject PRJNA573931, the genome accessions in `resources/references.tsv`, the relevant software publications and the associated manuscript when available. Record the repository commit used for an analysis; no manuscript DOI is assigned here.

Code is distributed under the MIT License. Public data remain subject to their source repositories' terms.
