# EVOKE Step 0: combine the three Pt-10 test samples into one raw-count Seurat object.
# Run this script from the repository root.

library(Seurat)

sample_id <- c("Pt10_Ad", "Pt10_P", "Pt10_Li")
sample_type <- c("Health", "Primary", "Metastasis")
GSM <- c("GSM7669805", "GSM7669807", "GSM7669806")
tissue <- c("Adjacent normal mucosa", "Primary gastric tumor", "Liver metastasis")

temp_dir <- tempfile()
dir.create(temp_dir)
unzip("example_data.zip", exdir = temp_dir)

sample_objects <- list()
for (i in seq_along(sample_id)) {
  counts <- Read10X(file.path(temp_dir, "example_data", sample_id[i]))
  sample_objects[[i]] <- CreateSeuratObject(
    counts = counts,
    project = sample_id[i],
    min.cells = 0,
    min.features = 0
  )
  sample_objects[[i]]$sample_id <- sample_id[i]
  sample_objects[[i]]$patient <- "Pt-10"
  sample_objects[[i]]$patient_id <- "Pt-10"
  sample_objects[[i]]$cancer_type <- "STAD"
  sample_objects[[i]]$species <- "Homo sapiens"
  sample_objects[[i]]$sample_type <- sample_type[i]
  sample_objects[[i]]$GSM <- GSM[i]
  sample_objects[[i]]$tissue <- tissue[i]
}

Pt10 <- merge(sample_objects[[1]], y = sample_objects[-1])
Pt10$orig.ident <- Pt10$sample_id
dir.create("data", showWarnings = FALSE)
saveRDS(Pt10, "data/Pt10_merged_raw.rds")

message("Saved one combined Seurat object: data/Pt10_merged_raw.rds")
