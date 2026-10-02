# fthapsinum-bmr-chemecology

Pathogen-side RNA-seq analysis of *Fusarium thapsinum* during sorghum stalk colonisation.

This workflow analyses host lignin-genotype and water-treatment associations with fungal RNA abundance and relative gene expression. It includes competitive mapping, QuantSeq-aware gene counting, gene-set analyses, protein annotation, promoter comparisons and a fixed sample-exclusion sensitivity analysis.

## Source data

The source experiment is NCBI BioProject **PRJNA573931**, described by Khasin et al. (2021), DOI **10.1186/s12870-021-03149-5**. It includes Tx430 wild type, CAD-deficient `bmr6` and COMT-deficient `bmr12` sorghum under two watering treatments.

| Analysis set | Libraries | Purpose |
| --- | ---: | --- |
| `fusarium_3dai` | 24 | Inoculated three-genotype by two-water design |
| `background_3dai` | 21 | Run-matched PDB controls |
| `mapping_3dai` | 45 | Primary competitive mapping and counting |
| `context_13dai` | 29 | Fifteen inoculated and fourteen control libraries |
| `host_3dai` | 48 | Balanced host-eigengene comparison |
| `macrophomina_3dai` | 24 | Pathogen-signal feasibility analysis |
| `macrophomina_mapping_3dai` | 45 | Macrophomina libraries plus shared controls |

`resources/samples.tsv` records 129 library accessions, sequencing URLs, checksums and experimental factors. Twenty-eight libraries are retained for provenance without assignment to these downstream sets. Three published bmr6 PDB names omit a day token; their inferred 3-day assignment is labelled and is used only in the balanced host analysis, not the primary background panel. Genome accessions, analysis membership and functional terms are in `resources/`.

## Run

From the repository root:

```bash
mamba env create -f environment.yml
mamba activate fthapsinum-bmr-chemecology
snakemake complete_study --dry-run --cores 1
snakemake complete_study --cores 12
```

`complete_study` includes the main fungal analysis, the Macrophomina feasibility arm and the 22-versus-21-library sensitivity fits. Individual targets are available:

```bash
snakemake --cores 12                         # main analysis
snakemake macrophomina --cores 12            # optional pathogen feasibility
snakemake sample_sensitivity --cores 12      # fixed sample-exclusion refits
```

Paths are relative to the repository. `config.yaml` controls reference roles, output locations, thresholds and thread counts. The full workflow may download and process sequencing data. Do not start a full run merely to inspect an existing result archive.

## Outputs

| Directory | Content |
| --- | --- |
| `results/fusarium/mapping_3dai/` | Competitive assignments and explicit pre/post-filter RNA metrics |
| `results/fusarium/count_model/` | QuantSeq counting-model selection |
| `results/fusarium/gene_qc/` | Library depth and background-locus filtering |
| `results/fusarium/expression/` | Gene-level and competitive gene-set results |
| `results/fusarium/module_audit/` | Fixed-set scores, matched sets and omission diagnostics |
| `results/fusarium/sample_sensitivity/` | Twenty-two versus twenty-one sample comparisons |
| `results/fusarium/sample_sensitivity_models/` | Four complete expression refits |
| `results/fusarium/focus_annotation/` | Protein-annotation evidence |
| `results/fusarium/regulatory_architecture/` | Locus positions, coexpression and complete-family motif tests |
| `results/fusarium/motif_conservation/` | Labelled comparative promoter follow-up |
| `results/fusarium/context_13dai/` | Later-sampling context |
| `results/host_module_coupling/` | Host-eigengene associations |
| `results/macrophomina/feasibility/` | Pathogen-signal feasibility |

Outputs are numerical tables, sequence files and model/session summaries. Raw reads, assemblies, alignments and generated results are not tracked. No figure-generation code or manuscript files are included.

## Analysis definitions

Both input-read and post-filter RPM are emitted for all arms. The historical reporting basis is explicitly recorded; the tissue-output offset remains post-filter primary reads. The 22-library primary analysis is preserved, with the named low-specificity library excluded only in the 21-library sensitivity fits. The well-watered bmr6-versus-wild-type difference persists in those fits, whereas interaction support weakens. Details and exact commands are in [RNA metrics and sample sensitivity](docs/rna_metrics_and_sensitivity.md).

Promoter discovery jointly adjusts all 10,272 reverse-complement-collapsed A/C/G/T 6- and 7-mers, without foreground-frequency filtering. No motif passes the configured full-family threshold. ATGCTAA/TTAGCAT has full-family q = 0.266708 and is followed only as an explicitly configured exploratory sequence. `auto` mode leaves zero or multiple passing candidates unselected. [Promoter analysis](docs/promoter_analysis.md) describes the testing family, selection modes and the distinction between measurable same-scaffold distances and unknown cross-scaffold distances.

The primary biological evidence is a within-experiment gene-set association. Neither reference assignment nor promoter conservation establishes co-infection, an inducing compound, transcription-factor binding or a virulence mechanism.

## Tests

```bash
python -m compileall -q workflow/scripts tests
pytest -q
snakemake complete_study --dry-run --cores 1
```

GitHub Actions checks syntax, unit tests and the main, optional pathogen and complete-study workflow plans. Numerical refits and corrected promoter analyses have separately been exercised on the archived local study inputs. CI does not reprocess the raw sequencing experiment.

## Citation and licence

Software authors are listed in `CITATION.cff`. Cite the source experiment, its BioProject, genome accessions, relevant methods and the exact repository version used. No manuscript DOI is assigned here. Code is distributed under the MIT License; public datasets retain their source terms.
