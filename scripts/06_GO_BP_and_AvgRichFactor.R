# EVOKE Step 6: GO-BP enrichment and AvgRichFactor.
# Run after Step 4.1 from the repository root.

library(clusterProfiler)
source("config/settings.R")
dir.create("outputs/patient_cellfunc", recursive = TRUE, showWarnings = FALSE)

gene_files <- list.files("outputs/patient_key_gene",
                         pattern = "_key_gene\\.tsv$",
                         full.names = TRUE)
patient_go <- data.frame(
  ID = character(), Description = character(), GeneRatio = character(),
  BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
  qvalue = numeric(), geneID = character(), Count = integer(),
  RichFactor = numeric(), significant = logical(), patient = character(),
  cancer_type = character(), celltype = character(), stage = character()
)

for (i in gene_files) {
  patientid <- sub("_key_gene\\.tsv$", "", basename(i))
  genes <- read.delim(i, stringsAsFactors = FALSE)
  genes <- genes[genes$event_category == "StageSpecific", ]
  if (nrow(genes) == 0) next

  SO <- readRDS(paste0("outputs/stage1/", patientid,
                       "_after_cell_annotation.rds"))
  if (SO$species[1] == "Homo sapiens") {
    OrgDb <- org.Hs.eg.db::org.Hs.eg.db
  } else {
    OrgDb <- org.Mm.eg.db::org.Mm.eg.db
  }
  groups <- unique(genes[, c("celltype", "stage")])

  for (j in seq_len(nrow(groups))) {
    selected <- genes$celltype == groups$celltype[j] &
      genes$stage == groups$stage[j]
    gene_list <- unique(genes$gene[selected])
    go <- enrichGO(
      gene = gene_list,
      OrgDb = OrgDb,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    )
    go <- as.data.frame(go)
    if (nrow(go) == 0) next

    go$RichFactor <- go$Count / as.numeric(sub(".*/", "", go$BgRatio))
    go$significant <- go$p.adjust < GO_BP_ADJ_P_CUTOFF
    go$patient <- SO$patient[1]
    go$cancer_type <- SO$cancer_type[1]
    go$celltype <- groups$celltype[j]
    go$stage <- groups$stage[j]
    patient_go <- rbind(patient_go, go)
  }
}

write.table(patient_go,
            "outputs/patient_cellfunc/patient_stage_specific_GO_BP.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

significant <- patient_go[patient_go$significant %in% TRUE, ]
avg_rich <- data.frame(
  cancer_type = character(), celltype = character(), stage = character(),
  ID = character(), Description = character(), n_patient_results = integer(),
  AvgRichFactor = numeric()
)
if (nrow(significant) > 0) {
  groups <- unique(significant[, c("cancer_type", "celltype", "stage",
                                    "ID", "Description")])
  for (i in seq_len(nrow(groups))) {
    selected <- significant$cancer_type == groups$cancer_type[i] &
      significant$celltype == groups$celltype[i] &
      significant$stage == groups$stage[i] &
      significant$ID == groups$ID[i]
    rows <- significant[selected, ]
    avg_rich <- rbind(avg_rich, data.frame(
      cancer_type = groups$cancer_type[i],
      celltype = groups$celltype[i],
      stage = groups$stage[i],
      ID = groups$ID[i],
      Description = groups$Description[i],
      n_patient_results = length(unique(rows$patient)),
      AvgRichFactor = mean(rows$RichFactor)
    ))
  }
}

write.table(avg_rich,
            "outputs/patient_cellfunc/AvgRichFactor_cancer_celltype_stage_GO.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Step 6 complete.")
