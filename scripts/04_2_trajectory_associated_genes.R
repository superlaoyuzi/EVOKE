# EVOKE Step 4.2: tradeSeq analysis of genes changing along each lineage.
# Run after Steps 1 and 3 from the repository root.


library(Seurat)
library(slingshot)
library(tradeSeq)
source("config/settings.R")
set.seed(RANDOM_SEED)

dir.create("outputs/sceGAM", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/trajectory_gene", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/trajectory_gene_umap", recursive = TRUE, showWarnings = FALSE)

file <- list.files("outputs/stage1",
                   pattern = "_after_cell_annotation\\.rds$",
                   full.names = TRUE)

for (file_index in seq_along(file)) {
  patientid <- sub("_after_cell_annotation\\.rds$", "",
                   basename(file[file_index]))
  message("[", file_index, "/", length(file), "] Processing patient: ", patientid)

  SO <- readRDS(file[file_index])
  curves <- readRDS(paste0("outputs/patient_slingshot/",
                           patientid, "_curves.rds"))
  BPPARAM <- BiocParallel::SnowParam(
    workers = TRADESEQ_WORKERS,
    type = "SOCK",
    RNGseed = RANDOM_SEED + file_index,
    progressbar = TRUE
  )

  counts <- GetAssayData(SO, assay = "RNA", layer = "counts")
  expressed <- Matrix::rowSums(counts > 0) / ncol(counts) >=
    TRADESEQ_MIN_CELL_FRACTION
  counts <- counts[expressed, , drop = FALSE]

  pseudotime <- slingPseudotime(curves, na = FALSE)
  cellWeights <- slingCurveWeights(curves)
  pseudotime <- pseudotime[colnames(counts), , drop = FALSE]
  cellWeights <- cellWeights[colnames(counts), , drop = FALSE]

  knot_file <- paste0("outputs/sceGAM/", patientid, "_evaluateK.csv")
  if (TRADESEQ_RUN_KNOT_CHECK && !file.exists(knot_file)) {
    # This diagnostic does not change the fixed six-knot model below.
    knot_check <- evaluateK(
      counts = counts,
      pseudotime = pseudotime,
      cellWeights = cellWeights,
      k = 3:10,
      nGenes = min(200, nrow(counts)),
      plot = FALSE,
      verbose = FALSE,
      parallel = TRUE,
      BPPARAM = BPPARAM
    )
    write.csv(knot_check, knot_file)
  }

  fit_file <- paste0("outputs/sceGAM/", patientid, "_GAM.rds")
  if (file.exists(fit_file)) {
    # Delete this cached file if the input object or tradeSeq settings changed.
    message("Using saved tradeSeq fit for ", patientid)
    fit <- readRDS(fit_file)
  } else {
    fit <- fitGAM(
      counts = counts,
      pseudotime = pseudotime,
      cellWeights = cellWeights,
      nknots = TRADESEQ_NKNOTS,
      verbose = FALSE,
      sce = TRUE,
      parallel = TRUE,
      BPPARAM = BPPARAM
    )
    saveRDS(fit, fit_file)
  }

  association <- associationTest(fit, global = FALSE, lineages = TRUE)
  start_end <- startVsEndTest(fit, global = FALSE, lineages = TRUE)
  expression <- GetAssayData(SO, assay = "RNA", layer = "data")
  umap <- Embeddings(SO, "umap")

  for (j in seq_len(ncol(pseudotime))) {
    lineage <- colnames(pseudotime)[j]
    association_p <- association[[paste0("pvalue_", j)]]
    start_end_p <- start_end[[paste0("pvalue_lineage", j)]]
    log2FC <- start_end[[paste0("logFClineage", j)]]

    result <- data.frame(
      gene = rownames(association),
      association_wald = association[[paste0("waldStat_", j)]],
      association_p = association_p,
      association_FDR = p.adjust(association_p, method = "BH"),
      start_end_wald = start_end[[paste0("waldStat_lineage", j)]],
      start_end_p = start_end_p,
      start_end_FDR = p.adjust(start_end_p, method = "BH"),
      log2FC = log2FC
    )
    result$trajectory_associated <-
      !is.na(result$association_FDR) &
      !is.na(result$start_end_FDR) &
      result$association_FDR < TRADESEQ_P_CUTOFF &
      result$start_end_FDR < TRADESEQ_P_CUTOFF &
      abs(result$log2FC) > TRADESEQ_LOG2FC_CUTOFF
    result$trajectory <- lineage
    result$patient <- SO$patient[1]
    result$cancer_type <- SO$cancer_type[1]
    result$test_method <- "tradeSeq Wald tests"

    write.table(result,
                paste0("outputs/trajectory_gene/", patientid, "_",
                       lineage, "_trajectory_gene.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    selected_genes <- result$gene[result$trajectory_associated]
    if (length(selected_genes) > 0) {
      lineage_cells <- rownames(cellWeights)[cellWeights[, j] > 0]
      if (length(lineage_cells) > 0) {
        gene_expression <- as.data.frame(
          t(as.matrix(expression[selected_genes, lineage_cells, drop = FALSE]))
        )
        gene_expression$barcode <- rownames(gene_expression)
        gene_expression$UMAP_1 <- umap[gene_expression$barcode, 1]
        gene_expression$UMAP_2 <- umap[gene_expression$barcode, 2]
        gene_expression$pseudotime <- pseudotime[gene_expression$barcode, j]
        gene_expression$cell_weight <- cellWeights[gene_expression$barcode, j]
        write.csv(gene_expression,
                  paste0("outputs/trajectory_gene_umap/", patientid, "_",
                         lineage, "_trajectory_gene_expression.csv"),
                  row.names = FALSE)
      }
    }
  }
}

message("Step 4.2 complete.")
