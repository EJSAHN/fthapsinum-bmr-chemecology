#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: fit_context_models.R COUNTS SAMPLES GENE_FILTER FOCUS_GENES OUTPUT_DIR MIN_COUNT")
}
counts_path <- args[[1]]
samples_path <- args[[2]]
filter_path <- args[[3]]
focus_path <- args[[4]]
out_root <- args[[5]]
minimum_count <- as.numeric(args[[6]])
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path, ...) read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, ...)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

counts <- as.matrix(read_tsv(counts_path, row.names = 1))
storage.mode(counts) <- "integer"
samples <- read_tsv(samples_path)
filter_table <- read_tsv(filter_path)
focus_table <- read_tsv(focus_path)
if (!all(samples$sample_id %in% colnames(counts))) stop("Sample/count mismatch")
counts <- counts[, samples$sample_id, drop = FALSE]
samples$group <- factor(
  paste(samples$genotype, samples$water, sep = "_"),
  levels = c("wt_dry", "wt_wet", "bmr12_dry", "bmr12_wet")
)
if (any(is.na(samples$group))) stop("Unexpected genotype-water group")

usable <- filter_table$gene_id[tolower(as.character(filter_table$usable_primary)) %in% c("true", "1", "yes")]
focus_ids <- intersect(focus_table$gene_id, rownames(counts))
keep <- filterByExpr(DGEList(counts = counts), group = samples$group, min.count = minimum_count, min.total.count = 10)
keep <- rownames(counts) %in% usable & keep
keep <- keep | (rownames(counts) %in% focus_ids & rowSums(counts) >= 10)
counts <- counts[keep, , drop = FALSE]
if (length(intersect(focus_ids, rownames(counts))) < 5) stop("Fewer than five focus genes remain")

design <- model.matrix(~0 + group, data = samples)
colnames(design) <- levels(samples$group)
rownames(design) <- samples$sample_id
if (qr(design)$rank != ncol(design)) stop("Rank-deficient context design")
contrasts <- makeContrasts(
  bmr12_vs_wt_dry = bmr12_dry - wt_dry,
  bmr12_vs_wt_wet = bmr12_wet - wt_wet,
  wet_vs_dry_wt = wt_wet - wt_dry,
  wet_vs_dry_bmr12 = bmr12_wet - bmr12_dry,
  bmr12_water_interaction = (bmr12_wet - bmr12_dry) - (wt_wet - wt_dry),
  levels = design
)

run_model <- function(name, scale) {
  output_dir <- file.path(out_root, name)
  dir.create(file.path(output_dir, "de"), recursive = TRUE, showWarnings = FALSE)
  if (scale == "state") {
    y <- DGEList(counts = counts)
    y <- calcNormFactors(y, method = "TMMwsp")
  } else {
    library_sizes <- as.numeric(samples$post_fastp_reads)
    if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) stop("Invalid tissue read counts")
    y <- DGEList(counts = counts, lib.size = library_sizes)
    y$samples$norm.factors <- 1
  }
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  logcpm <- cpm(y, log = TRUE, prior.count = 2)
  write_tsv(data.frame(gene_id = rownames(logcpm), logcpm, check.names = FALSE), file.path(output_dir, "logCPM.tsv"))

  summaries <- list()
  for (contrast_name in colnames(contrasts)) {
    test <- glmQLFTest(fit, contrast = contrasts[, contrast_name])
    table <- topTags(test, n = Inf, sort.by = "none")$table
    table$gene_id <- rownames(table)
    rownames(table) <- NULL
    table <- table[, c("gene_id", setdiff(colnames(table), "gene_id"))]
    write_tsv(table, file.path(output_dir, "de", paste0(safe_name(contrast_name), ".tsv")))
    summaries[[contrast_name]] <- data.frame(
      model = name,
      scale = scale,
      contrast = contrast_name,
      genes = nrow(table),
      fdr05 = sum(table$FDR < 0.05, na.rm = TRUE),
      min_fdr = min(table$FDR, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  write_tsv(do.call(rbind, summaries), file.path(output_dir, "contrast_summary.tsv"))

  indices <- ids2indices(
    list(aromatic_phenolic_metabolism = intersect(focus_ids, rownames(logcpm))),
    rownames(logcpm), remove.empty = TRUE
  )
  v <- voom(y, design = design, plot = FALSE, normalize.method = "none")
  camera_rows <- list()
  for (contrast_name in colnames(contrasts)) {
    result <- camera(v$E, index = indices, design = design, contrast = contrasts[, contrast_name], sort = FALSE)
    result$module <- rownames(result)
    result$contrast <- contrast_name
    result$model <- name
    rownames(result) <- NULL
    camera_rows[[contrast_name]] <- result
  }
  write_tsv(do.call(rbind, camera_rows), file.path(output_dir, "camera.tsv"))
  write_json(
    list(model = name, scale = scale, samples = ncol(y), genes = nrow(y), groups = as.list(table(samples$group))),
    file.path(output_dir, "model.json"), pretty = TRUE, auto_unbox = TRUE
  )
}

run_model("state", "state")
run_model("tissue", "tissue")
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
