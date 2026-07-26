#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop("Usage: fit_expression_models.R COUNTS_PRIMARY SAMPLES_PRIMARY COUNTS_ALL SAMPLES_ALL GENE_SETS MAPPING OUTPUT_DIR FDR MIN_ABS_LOG2FC")
}
counts_primary_path <- args[[1]]
samples_primary_path <- args[[2]]
counts_all_path <- args[[3]]
samples_all_path <- args[[4]]
gene_sets_path <- args[[5]]
mapping_path <- args[[6]]
out_root <- args[[7]]
fdr_cutoff <- as.numeric(args[[8]])
lfc_cutoff <- as.numeric(args[[9]])
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

read_counts <- function(path) {
  x <- read.delim(path, check.names = FALSE, row.names = 1)
  as.matrix(x)
}
read_samples <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  x$group <- factor(
    paste(x$genotype, x$water, sep = "_"),
    levels = c("wt_dry", "bmr6_dry", "bmr12_dry", "wt_wet", "bmr6_wet", "bmr12_wet")
  )
  x
}
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

counts_primary <- read_counts(counts_primary_path)
counts_all <- read_counts(counts_all_path)
samples_primary <- read_samples(samples_primary_path)
samples_all <- read_samples(samples_all_path)
mapping <- read.delim(mapping_path, check.names = FALSE, stringsAsFactors = FALSE)
rownames(mapping) <- mapping$sample_id

sets_table <- read.delim(gene_sets_path, check.names = FALSE, stringsAsFactors = FALSE)
functional_sets <- split(sets_table$gene_id, sets_table$module)

contrast_matrix <- function(columns) {
  makeContrasts(
    bmr6_vs_wt_dry = bmr6_dry - wt_dry,
    bmr12_vs_wt_dry = bmr12_dry - wt_dry,
    bmr6_vs_wt_wet = bmr6_wet - wt_wet,
    bmr12_vs_wt_wet = bmr12_wet - wt_wet,
    wet_vs_dry_wt = wt_wet - wt_dry,
    wet_vs_dry_bmr6 = bmr6_wet - bmr6_dry,
    wet_vs_dry_bmr12 = bmr12_wet - bmr12_dry,
    bmr6_water_interaction = (bmr6_wet - bmr6_dry) - (wt_wet - wt_dry),
    bmr12_water_interaction = (bmr12_wet - bmr12_dry) - (wt_wet - wt_dry),
    bmr12_vs_bmr6_dry = bmr12_dry - bmr6_dry,
    bmr12_vs_bmr6_wet = bmr12_wet - bmr6_wet,
    levels = columns
  )
}

make_y <- function(counts, samples, scale) {
  counts <- counts[, samples$sample_id, drop = FALSE]
  if (scale == "state") {
    y <- DGEList(counts = counts)
    y <- calcNormFactors(y, method = "TMMwsp")
  } else {
    tissue_reads <- mapping[samples$sample_id, "post_fastp_reads"]
    if (any(!is.finite(tissue_reads)) || any(tissue_reads <= 0)) stop("Invalid tissue read counts")
    y <- DGEList(counts = counts, lib.size = tissue_reads)
    y$samples$norm.factors <- 1
  }
  y
}

run_camera <- function(y, design, contrasts, output_dir) {
  v <- voom(y, design = design, plot = FALSE, normalize.method = "none")
  sets <- lapply(functional_sets, function(ids) intersect(ids, rownames(v$E)))
  sets <- sets[lengths(sets) >= 5 & lengths(sets) <= 500]
  if (length(sets) == 0) return(data.frame())
  indices <- ids2indices(sets, rownames(v$E), remove.empty = TRUE)
  rows <- list()
  for (contrast in colnames(contrasts)) {
    result <- camera(v$E, index = indices, design = design, contrast = contrasts[, contrast], sort = FALSE)
    result$module <- rownames(result)
    result$contrast <- contrast
    rownames(result) <- NULL
    rows[[contrast]] <- result
  }
  combined <- do.call(rbind, rows)
  rownames(combined) <- NULL
  write_tsv(combined, file.path(output_dir, "camera.tsv"))
  combined
}

run_model <- function(name, counts, samples, scale) {
  output_dir <- file.path(out_root, name)
  dir.create(file.path(output_dir, "de"), recursive = TRUE, showWarnings = FALSE)
  y <- make_y(counts, samples, scale)
  design <- model.matrix(~0 + group, data = samples)
  colnames(design) <- sub("^group", "", colnames(design))
  rownames(design) <- samples$sample_id
  if (qr(design)$rank != ncol(design)) stop(paste("Rank-deficient design:", name))
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  contrasts <- contrast_matrix(colnames(design))
  logcpm <- cpm(y, log = TRUE, prior.count = 2)
  write_tsv(data.frame(gene_id = rownames(logcpm), logcpm, check.names = FALSE), file.path(output_dir, "logCPM.tsv"))
  write_tsv(data.frame(sample_id = rownames(design), design, check.names = FALSE), file.path(output_dir, "design.tsv"))

  summaries <- list()
  for (contrast in colnames(contrasts)) {
    test <- glmQLFTest(fit, contrast = contrasts[, contrast])
    table <- topTags(test, n = Inf, sort.by = "none")$table
    table$gene_id <- rownames(table)
    rownames(table) <- NULL
    table$significant <- table$FDR < fdr_cutoff & abs(table$logFC) >= lfc_cutoff
    table <- table[, c("gene_id", setdiff(colnames(table), "gene_id"))]
    write_tsv(table, file.path(output_dir, "de", paste0(safe_name(contrast), ".tsv")))
    summaries[[contrast]] <- data.frame(
      model = name, scale = scale, contrast = contrast, genes = nrow(table),
      fdr05 = sum(table$FDR < fdr_cutoff, na.rm = TRUE),
      fdr05_lfc = sum(table$significant, na.rm = TRUE),
      min_fdr = min(table$FDR, na.rm = TRUE), stringsAsFactors = FALSE
    )
  }
  summary <- do.call(rbind, summaries)
  write_tsv(summary, file.path(output_dir, "contrast_summary.tsv"))

  design_global <- model.matrix(~group, data = samples)
  y_global <- estimateDisp(y, design_global, robust = TRUE)
  fit_global <- glmQLFit(y_global, design_global, robust = TRUE)
  global_test <- glmQLFTest(fit_global, coef = 2:ncol(design_global))
  global_table <- topTags(global_test, n = Inf, sort.by = "none")$table
  global_table$gene_id <- rownames(global_table)
  rownames(global_table) <- NULL
  write_tsv(global_table[, c("gene_id", setdiff(colnames(global_table), "gene_id"))], file.path(output_dir, "global_group_test.tsv"))

  design_interaction <- model.matrix(~genotype * water, data = samples)
  y_interaction <- estimateDisp(y, design_interaction, robust = TRUE)
  fit_interaction <- glmQLFit(y_interaction, design_interaction, robust = TRUE)
  interaction_columns <- grep(":", colnames(design_interaction))
  if (length(interaction_columns) > 0) {
    interaction_test <- glmQLFTest(fit_interaction, coef = interaction_columns)
    interaction_table <- topTags(interaction_test, n = Inf, sort.by = "none")$table
    interaction_table$gene_id <- rownames(interaction_table)
    rownames(interaction_table) <- NULL
    write_tsv(interaction_table[, c("gene_id", setdiff(colnames(interaction_table), "gene_id"))], file.path(output_dir, "global_interaction_test.tsv"))
  }

  camera_result <- run_camera(y, design, contrasts, output_dir)
  info <- list(
    name = name, scale = scale, samples = ncol(y), genes = nrow(y),
    common_dispersion = unname(y$common.dispersion),
    edgeR_version = as.character(packageVersion("edgeR")),
    limma_version = as.character(packageVersion("limma")),
    camera_sets = if (nrow(camera_result)) length(unique(camera_result$module)) else 0
  )
  write_json(info, file.path(output_dir, "model.json"), pretty = TRUE, auto_unbox = TRUE)
  summary
}

summaries <- list(
  run_model("primary_state", counts_primary, samples_primary, "state"),
  run_model("primary_tissue", counts_primary, samples_primary, "tissue"),
  run_model("all_state", counts_primary, samples_all, "state"),
  run_model("primary_state_inclusive", counts_all, samples_primary, "state")
)
combined <- do.call(rbind, summaries)
write_tsv(combined, file.path(out_root, "contrast_summary.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
