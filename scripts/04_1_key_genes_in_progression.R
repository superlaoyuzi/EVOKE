# EVOKE Step 4.1: identify stage-specific and cross-stage progression genes.
# Run after Step 1 from the repository root.

library(Seurat)
source("config/settings.R")

dir.create("outputs/patient_key_gene", recursive = TRUE, showWarnings = FALSE)
file <- list.files("outputs/stage1",
                   pattern = "_after_cell_annotation\\.rds$",
                   full.names = TRUE)
skipped_small_comparisons <- 0L

for (i in file) {
  SO <- readRDS(i)
  patientid <- sub("_after_cell_annotation\\.rds$", "", basename(i))
  stages <- STAGE_ORDER[STAGE_ORDER %in% unique(as.character(SO$sample_type))]
  cell_types <- unique(as.character(SO$celltype))

  events <- data.frame(
    event_category = character(), gene = character(), celltype = character(),
    stage = character(), from_stage = character(), to_stage = character(),
    log2FC = numeric(), p_value = numeric(), adjusted_p = numeric(),
    patient = character(), cancer_type = character()
  )

  # Within each stage, compare one cell type with all other cell types.
  for (stage in stages) {
    stage_SO <- subset(SO, subset = sample_type == stage)
    if (length(unique(stage_SO$celltype)) < 2) next
    Idents(stage_SO) <- "celltype"

    for (cell_type in unique(as.character(stage_SO$celltype))) {
      if (cell_type == "Unknown") next
      n_cell_type <- sum(as.character(stage_SO$celltype) == cell_type)
      if (n_cell_type < 3L || nrow(stage_SO) - n_cell_type < 3L) {
        skipped_small_comparisons <- skipped_small_comparisons + 1L
        next
      }
      de <- FindMarkers(
        stage_SO, ident.1 = cell_type, test.use = "wilcox",
        min.pct = 0, logfc.threshold = 0,
        p.adjust.method = "bonferroni", verbose = FALSE
      )
      de <- de[!is.na(de$p_val_adj) & de$p_val_adj < 0.05 &
                 de$avg_log2FC > 0.5, , drop = FALSE]
      if (nrow(de) == 0) next
      events <- rbind(events, data.frame(
        event_category = "StageSpecific",
        gene = rownames(de),
        celltype = cell_type,
        stage = stage,
        from_stage = NA,
        to_stage = NA,
        log2FC = de$avg_log2FC,
        p_value = de$p_val,
        adjusted_p = de$p_val_adj,
        patient = SO$patient[1],
        cancer_type = SO$cancer_type[1]
      ))
    }
  }

  # Within each cell type, compare every pair of consecutive observed stages.
  # Keep a gene only when every comparison passes and its direction is consistent.
  for (cell_type in cell_types) {
    if (cell_type == "Unknown") next
    cell_SO <- subset(SO, subset = celltype == cell_type)
    present <- STAGE_ORDER[STAGE_ORDER %in% unique(as.character(cell_SO$sample_type))]
    if (length(present) < 3) next
    stage_cell_counts <- table(factor(as.character(cell_SO$sample_type),
                                      levels = present))
    if (any(stage_cell_counts < 3L)) {
      skipped_small_comparisons <- skipped_small_comparisons + length(present) - 1L
      next
    }
    Idents(cell_SO) <- "sample_type"

    de_list <- list()
    for (k in seq_len(length(present) - 1)) {
      de_list[[k]] <- FindMarkers(
        cell_SO, ident.1 = present[k + 1], ident.2 = present[k],
        test.use = "wilcox", min.pct = 0, logfc.threshold = 0,
        p.adjust.method = "bonferroni", verbose = FALSE
      )
    }
    common_genes <- Reduce(intersect, lapply(de_list, rownames))

    for (gene in common_genes) {
      log2FC <- numeric(length(de_list))
      adjusted_p <- numeric(length(de_list))
      p_value <- numeric(length(de_list))
      for (k in seq_along(de_list)) {
        log2FC[k] <- de_list[[k]][gene, "avg_log2FC"]
        adjusted_p[k] <- de_list[[k]][gene, "p_val_adj"]
        p_value[k] <- de_list[[k]][gene, "p_val"]
      }
      if (!anyNA(adjusted_p) && !anyNA(log2FC) &&
          all(adjusted_p < 0.05) && all(abs(log2FC) > 0.5) &&
          (all(log2FC > 0) || all(log2FC < 0))) {
        for (k in seq_along(de_list)) {
          events <- rbind(events, data.frame(
            event_category = "CrossStageConsistent",
            gene = gene,
            celltype = cell_type,
            stage = NA,
            from_stage = present[k],
            to_stage = present[k + 1],
            log2FC = log2FC[k],
            p_value = p_value[k],
            adjusted_p = adjusted_p[k],
            patient = SO$patient[1],
            cancer_type = SO$cancer_type[1]
          ))
        }
      }
    }
  }

  write.table(events,
              paste0("outputs/patient_key_gene/", patientid, "_key_gene.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

message("Step 4.1 complete. Comparisons skipped because a group had fewer than 3 cells: ",
        skipped_small_comparisons)
