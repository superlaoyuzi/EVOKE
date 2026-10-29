# EVOKE Step 2: cell proportions, Ro/e tissue preference and Augur scores.
# Run after Step 1 from the repository root.

library(Seurat)
library(Augur)
source("config/settings.R")

dir.create("outputs/patient_cell_prop", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/patient_Roe", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/patient_augur", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/patient_celltable", recursive = TRUE, showWarnings = FALSE)

file <- list.files("outputs/stage1",
                   pattern = "_after_cell_annotation\\.rds$",
                   full.names = TRUE)

for (i in file) {
  SO <- readRDS(i)
  patientid <- sub("_after_cell_annotation\\.rds$", "", basename(i))

  cell_prop <- prop.table(table(SO$celltype, SO$sample_type), margin = 2)
  cell_prop <- as.data.frame(as.table(cell_prop))
  colnames(cell_prop) <- c("celltype", "sample_type", "proportion")
  write.csv(cell_prop,
            paste0("outputs/patient_cell_prop/", patientid, "_cell_prop.csv"),
            row.names = FALSE)

  observed <- table(SO$celltype, SO$sample_type)
  expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)
  Roe <- observed / expected
  Roe <- as.data.frame(as.table(Roe))
  colnames(Roe) <- c("celltype", "sample_type", "RoE")
  Roe$status <- ifelse(Roe$RoE > ROE_ENRICHMENT_CUTOFF, "enriched",
                       ifelse(Roe$RoE < ROE_DEPLETION_CUTOFF,
                              "depleted", "neither"))
  write.csv(Roe,
            paste0("outputs/patient_Roe/", patientid, "_Roe.csv"),
            row.names = FALSE)

  # Augur defaults: 50 balanced subsamples and 20 cells per state.
  augur <- calculate_auc(SO, cell_type_col = "celltype",
                         label_col = "sample_type", n_threads = 8)
  augur$AUC$high_confidence <-
    augur$AUC$auc > AUGUR_HIGH_CONFIDENCE_CUTOFF
  auc <- as.data.frame(augur$AUC)
  names(auc)[names(auc) == "cell_type"] <- "celltype"
  names(auc)[names(auc) == "auc"] <- "perturbation_score"
  write.csv(augur$AUC,
            paste0("outputs/patient_augur/", patientid, "_augur_auc.csv"),
            row.names = FALSE)
  saveRDS(augur,
          paste0("outputs/patient_augur/", patientid, "_augur.rds"))

  cell_table <- merge(cell_prop, Roe, by = c("celltype", "sample_type"),
                      all = TRUE)
  cell_table <- merge(cell_table, auc, by = "celltype", all = TRUE)
  cell_table$patient <- SO$patient[1]
  cell_table$cancer_type <- SO$cancer_type[1]
  write.csv(cell_table,
            paste0("outputs/patient_celltable/", patientid, "_celltable.csv"),
            row.names = FALSE)
}

message("Step 2 complete.")
