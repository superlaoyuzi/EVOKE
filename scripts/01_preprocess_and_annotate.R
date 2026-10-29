# EVOKE Step 1: QC, sample-wise normalization, clustering, UMAP and scType annotation.
# Run after scripts/00_prepare_pt10.R from the repository root.

library(Seurat)
library(dplyr)
library(HGNChelper)
library(openxlsx)
library(scales)
source("config/settings.R")
set.seed(RANDOM_SEED)

Pt10 <- readRDS("data/Pt10_merged_raw.rds")
sample_names <- unique(Pt10$sample_id)
sample_objects <- list()

for (i in seq_along(sample_names)) {
  x <- subset(Pt10, subset = sample_id == sample_names[i])
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^(MT-|mt-)")
  x <- subset(x, subset = nFeature_RNA >= MIN_FEATURES &
                percent.mt <= MAX_MITO_PERCENT)
  x <- NormalizeData(x, normalization.method = "LogNormalize",
                     scale.factor = NORMALIZATION_SCALE, verbose = FALSE)
  sample_objects[[i]] <- x
}

SO <- merge(sample_objects[[1]], y = sample_objects[-1], merge.data = TRUE)
SO <- JoinLayers(SO, assay = "RNA")
DefaultAssay(SO) <- "RNA"
SO <- FindVariableFeatures(SO, selection.method = "vst",
                           nfeatures = N_VARIABLE_FEATURES, verbose = FALSE)
SO <- ScaleData(SO, verbose = FALSE)
SO <- RunPCA(SO, features = VariableFeatures(SO), npcs = N_PCS,
             verbose = FALSE)
SO <- FindNeighbors(SO, dims = 1:N_PCS, verbose = FALSE)
SO <- FindClusters(SO, resolution = CLUSTER_RESOLUTION, verbose = FALSE)
SO$trajectory_cluster <- as.character(SO$seurat_clusters)
SO <- RunUMAP(SO, dims = 1:N_PCS, seed.use = RANDOM_SEED, verbose = FALSE)

source(SCTYPE_WRAPPER)
SO <- run_sctype(SO, known_tissue_type = SCTYPE_TISSUE,
                 custom_marker_file = SCTYPE_MARKERS,
                 name = "sctype_classification", plot = FALSE)
SO <- run_sctype(SO, known_tissue_type = "Immune system",
                 custom_marker_file = SCTYPE_MARKERS,
                 name = "sctype_immune_classification", plot = FALSE)

SO$celltype <- as.character(SO$sctype_classification)
immune <- SO$celltype %in% c("Lymphoid cells", "Myeloid cells", "Unknown")
SO$celltype[immune] <- as.character(SO$sctype_immune_classification[immune])

dir.create("outputs/stage1", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/patient_umap", recursive = TRUE, showWarnings = FALSE)
review_file <- "outputs/stage1/STAD_Pt-10_cluster_celltype_review.csv"
if (!file.exists(review_file)) {
  review <- unique(data.frame(
    trajectory_cluster = SO$trajectory_cluster,
    automatic_label = SO$celltype
  ))
  review$reviewed_label <- ""
  write.csv(review, review_file, row.names = FALSE)
}
review <- read.csv(review_file, stringsAsFactors = FALSE)
reviewed_label <- review$reviewed_label[
  match(SO$trajectory_cluster, as.character(review$trajectory_cluster))
]
corrected <- !is.na(reviewed_label) & nzchar(trimws(reviewed_label))
SO$celltype[corrected] <- reviewed_label[corrected]
tumor_epithelial <- SO$sample_type %in% c("Primary", "Metastasis") &
  grepl("epithelial", SO$celltype, ignore.case = TRUE)
SO$celltype[tumor_epithelial] <- "Cancer cells"
umap <- as.data.frame(Embeddings(SO, "umap"))
umap$celltype <- SO$celltype
umap$sample_type <- SO$sample_type
write.csv(umap, "outputs/patient_umap/STAD_Pt-10_umap.csv")

saveRDS(SO, "outputs/stage1/STAD_Pt-10_after_cell_annotation.rds")
message("Step 1 complete: outputs/stage1/STAD_Pt-10_after_cell_annotation.rds")
