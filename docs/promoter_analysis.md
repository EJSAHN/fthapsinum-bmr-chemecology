# Promoter analysis

## Hypothesis family

The configured word lengths are 6 and 7. A word and its reverse complement define one hypothesis. The complete A/C/G/T family contains 2,080 canonical 6-mers and 8,192 canonical 7-mers. A motif contributes at most one presence per promoter, irrespective of its number of sites. Windows containing ambiguous bases are skipped.

For every word, a one-sided Fisher test compares promoter presence in the target and matched-background groups. Benjamini-Hochberg adjustment is applied jointly to all 10,272 hypotheses. Zero-target and unobserved words remain in that family. Counts are not used to select a smaller family before adjustment. `fisher_fdr` is a compatibility alias for `bh_full_family_q`, not a second correction.

`motif_search_summary.json` records the family definition, threshold and number of passing hypotheses. `promoter_kmer_enrichment.tsv` contains the full family, its size and each word's passing flag.

## Candidate selection

`auto` mode selects a candidate only when exactly one motif passes the complete-family threshold. Zero passing motifs produce `NO_DISCOVERY_MOTIF`; multiple passing motifs produce an unselected summary. Neither case silently selects the smallest P value. Unselected runs clear comparative table outputs so that a result from an earlier run cannot be mistaken for a current result.

`exploratory` mode requires an explicit motif. It returns the complete-family q value with `passes_full_family_fdr` and `regulatory_function_validated` fields. This mode permits descriptive study of a failed discovery candidate without changing its inferential status.

In the study inputs, ATGCTAA/TTAGCAT was present in six of ten target promoters and 38 of 557 matched-background promoters. Its odds ratio was 20.4868, one-sided Fisher P was 2.59645 x 10^-5 and full-family q was 0.266708. No word passed q <= 0.10. The supplied configuration therefore treats this sequence as exploratory.

Fixed-candidate matched-set tail probabilities ask about a specified sequence. They do not repeat the complete discovery-and-selection procedure within each random set and are not scan-wide adjusted probabilities. Conservation in homologous promoters is also not independent biological validation or evidence of transcription-factor binding.

## Projected locations

Location summaries retain sequence-record identifiers. Same-record distances are measurable; distances between separate records are unknown. The comparative layout table reports interval-to-interval gaps among same-record pairs and leaves cross-record gaps missing. The within-FL-4 clustering summaries separately use gene-centre distances, as labelled in that analysis.

The ten target loci occupy ten sequence records in each of the two F. thapsinum references. This does not determine their chromosome-level arrangement. Corresponding loci in related references provide descriptive comparative context, not a new orthology or synteny validation.

## Reproducibility scope

The affected production commands were tested against the archived study inputs. Complete-family statistics, regenerated promoter sequences, matching pools and comparative calls were compared with separately reviewed outputs. Protein-to-genome projections were reused. This targeted check did not rerun raw-read processing, RNA alignment, gene counting or edgeR models.
