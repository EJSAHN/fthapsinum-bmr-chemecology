# RNA metrics and sample-exclusion analysis

## RNA denominators

For every analysis arm, `summarize_mapping.py` reports two explicitly named metrics:

- `target_unique_rpm_input`: target-reference unique primary alignments per million reads before fastp filtering.
- `target_unique_rpm_post_filter`: the same numerator per million reads after filtering.

Both read counts come from the fastp JSON report. The compatibility column `target_unique_rpm` follows `--rpm-basis`, recorded in `rpm_basis`. The 3-day Fusarium comparison uses `input`; the 13-day Fusarium and Macrophomina contextual analyses retain `post_filter`. This preserves the denominators used in their original models. Use an explicitly named column when making a common-scale descriptive comparison between arms.

For Fusarium thapsinum, `target_specificity` is its fraction among reads uniquely assigned to it or the F. verticillioides, F. proliferatum and F. fujikuroi references. `target_fraction_all_fungal` additionally includes the other fungal reference groups. Comparator membership is recorded in `specificity_comparators`; missing counts are errors rather than implicit zeros. These are reference-assignment fractions, not species-identification or biomass assays.

The tissue-output expression offset remains the number of post-filter primary reads. Changing an RNA-burden display denominator does not change that offset.

## Fixed sample-exclusion comparison

`resources/sample_sensitivity.json` names one excluded library, `bmr6_fus_dry_3_d`, because of its low target-reference assignment fraction. The original 22-library primary analysis is retained. A separate 21-library fit uses the same 3,926 genes and the same functional-set membership. Neither genes nor excluded samples are selected by the resulting P values.

```bash
snakemake sample_sensitivity --cores 12
```

The target prepares fixed inputs, refits fungal-relative and tissue-output models for both subsets, and summarises the three focal score contrasts. The R command repeats TMMwsp normalisation, dispersion estimation and robust quasi-likelihood fitting. CAMERA uses the voom expression matrix with inter-gene correlation 0.01, matching the recorded analysis. It reports both within-contrast FDR and BH adjustment over all set-by-contrast pairs.

The Python summary recomputes gene-wise population-standardised scores from each refit. HC3 intervals and tests use the normal reference distribution. The water-label enumeration permutes water labels within genotype; it is not an interaction-only null that permits arbitrary water main effects. Fixed-score deletion diagnostics are kept separate from the full refits. Completion never depends on significance.

The locally executed study comparison retained the well-watered bmr6-versus-wild-type score difference (1.0653 versus 1.0632). After exclusion, the RNA-adjusted interaction had HC3 P = 0.05249 and within-contrast CAMERA FDR = 0.39266. Thus the pairwise host comparison and the interaction must be reported separately.

## Outputs

- `results/fusarium/sample_sensitivity/`: score contrasts, refitted scores, focal CAMERA comparisons, assignment QC and completion metadata.
- `results/fusarium/sample_sensitivity_models/`: four complete refits, including normalisation factors, logCPM, all eleven gene-level contrasts, CAMERA and package versions.
- `work/sample_sensitivity/inputs/`: fixed counts, sample design, gene universe and settings.

The raw-read workflow and these downstream refits are distinct execution scopes. The local sample-exclusion check did not re-align reads, recount genes or rerun promoter analysis.
