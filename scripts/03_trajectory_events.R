# EVOKE Step 3: Slingshot trajectories and progression-stage sliding windows.
# Run after Step 1 from the repository root.

library(Seurat)
library(slingshot)
library(zoo)
source("config/settings.R")
set.seed(RANDOM_SEED)

dir.create("outputs/patient_slingshot", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/patient_trajectory_table", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/lineage_slingshot_table", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/lineage_slingshot_plot", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/curve_umap_table", recursive = TRUE, showWarnings = FALSE)

file <- list.files("outputs/stage1",
                   pattern = "_after_cell_annotation\\.rds$",
                   full.names = TRUE)

for (i in file) {
  SO <- readRDS(i)
  patientid <- sub("_after_cell_annotation\\.rds$", "", basename(i))
  stages <- STAGE_ORDER[STAGE_ORDER %in% unique(as.character(SO$sample_type))]
  coordinates <- Embeddings(SO, "umap")

  cluster_stage <- table(SO$trajectory_cluster, SO$sample_type)
  start_stage_fraction <- cluster_stage[, stages[1]] / rowSums(cluster_stage)
  root_cluster <- names(which.max(start_stage_fraction))

  lineages <- as.SlingshotDataSet(getLineages(
    data = coordinates,
    clusterLabels = SO$trajectory_cluster,
    start.clus = root_cluster,
    dist.method = "slingshot"
  ))
  curves <- as.SlingshotDataSet(getCurves(
    lineages, thresh = 0.1, stretch = 2,
    allow.breaks = FALSE, approx_points = 150
  ))
  saveRDS(lineages,
          paste0("outputs/patient_slingshot/", patientid, "_slingshot_data.rds"))
  saveRDS(curves,
          paste0("outputs/patient_slingshot/", patientid, "_curves.rds"))
  write.csv(data.frame(
    patient = SO$patient[1],
    cancer_type = SO$cancer_type[1],
    start_stage = stages[1],
    root_cluster = root_cluster,
    start_stage_fraction = start_stage_fraction[root_cluster]
  ), paste0("outputs/patient_slingshot/", patientid, "_root.csv"),
  row.names = FALSE)

  pseudotime <- slingPseudotime(curves, na = TRUE)
  trajectory <- data.frame()

  for (j in seq_len(ncol(pseudotime))) {
    lineage <- colnames(pseudotime)[j]
    selected <- which(!is.na(pseudotime[, j]))
    selected <- selected[order(pseudotime[selected, j])]
    n_cells <- length(selected)
    window_size <- min(n_cells, ceiling(n_cells / 5000) * 100)
    sample_order <- as.character(SO$sample_type[selected])

    window_table <- data.frame(
      barcode = colnames(SO)[selected],
      pseudotime = pseudotime[selected, j],
      sample_type = sample_order
    )
    for (stage in stages) {
      window_table[[paste0("prop_", make.names(stage))]] <- zoo::rollapply(
        as.numeric(sample_order == stage), width = window_size,
        FUN = mean, fill = NA, align = "right"
      )
    }
    write.csv(window_table,
              paste0("outputs/lineage_slingshot_table/", patientid, "_",
                     lineage, "_frequency.csv"),
              row.names = FALSE)
    pdf(paste0("outputs/lineage_slingshot_plot/", patientid, "_",
               lineage, "_frequency.pdf"))
    matplot(seq_len(n_cells),
            as.matrix(window_table[, paste0("prop_", make.names(stages))]),
            type = "l", lty = 1, lwd = 2, ylim = c(0, 1),
            xlab = "Pseudotime order", ylab = "Stage frequency")
    legend("topright", legend = stages, col = seq_along(stages),
           lty = 1, lwd = 2, bty = "n")
    dev.off()

    curve_xy <- curves@curves[[j]]$s
    pdf(paste0("outputs/lineage_slingshot_plot/", patientid, "_",
               lineage, "_curve.pdf"))
    plot(coordinates[, 1], coordinates[, 2], pch = 16, cex = 0.4,
         col = "grey75", xlab = "UMAP 1", ylab = "UMAP 2")
    lines(curve_xy[, 1], curve_xy[, 2], lwd = 2)
    arrows(curve_xy[nrow(curve_xy) - 1, 1], curve_xy[nrow(curve_xy) - 1, 2],
           curve_xy[nrow(curve_xy), 1], curve_xy[nrow(curve_xy), 2],
           length = 0.08, lwd = 2)
    dev.off()

    write.csv(data.frame(
      barcode = colnames(SO),
      UMAP_1 = coordinates[, 1],
      UMAP_2 = coordinates[, 2],
      celltype = SO$celltype,
      sample_type = SO$sample_type,
      pseudotime = pseudotime[, j]
    ), paste0("outputs/curve_umap_table/", patientid, "_",
              lineage, "_curve_umap.csv"),
    row.names = FALSE)
    write.csv(curve_xy,
              paste0("outputs/curve_umap_table/", patientid, "_",
                     lineage, "_curve_line.csv"),
              row.names = FALSE)

    lineage_clusters <- slingLineages(curves)[[j]]
    trajectory <- rbind(trajectory, data.frame(
      trajectory = lineage,
      cell_trajectory = paste(lineage_clusters, collapse = "->"),
      cells_on_lineage = n_cells,
      sliding_window = window_size,
      StartStage = stages[1],
      EndStage = sample_order[n_cells],
      root_cluster = root_cluster,
      patient = SO$patient[1],
      cancer_type = SO$cancer_type[1]
    ))
  }

  write.csv(trajectory,
            paste0("outputs/patient_trajectory_table/", patientid,
                   "_trajectory_table.csv"),
            row.names = FALSE)
}

message("Step 3 complete.")
