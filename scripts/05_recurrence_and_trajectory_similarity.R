# EVOKE Step 5: recurrence scores and trajectory Jaccard similarity.
# Run after Steps 4.1 and 4.2 from the repository root.

source("config/settings.R")
dir.create("outputs/database_events", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/trajectory_similarity", recursive = TRUE, showWarnings = FALSE)

case_files <- list.files("outputs/stage1",
                         pattern = "_after_cell_annotation\\.rds$",
                         full.names = TRUE)
coverage <- data.frame()

for (i in case_files) {
  SO <- readRDS(i)
  case_id <- sub("_after_cell_annotation\\.rds$", "", basename(i))
  stages <- STAGE_ORDER[STAGE_ORDER %in% unique(as.character(SO$sample_type))]
  stage_table <- table(SO$celltype, SO$sample_type)
  coverage <- rbind(coverage, data.frame(
    case_id = case_id,
    patient = SO$patient[1],
    cancer_type = SO$cancer_type[1],
    StageSpecific = length(stages) >= 1,
    CrossStageConsistent = any(rowSums(stage_table > 0) >= 3)
  ))
}

event_files <- list.files("outputs/patient_key_gene",
                          pattern = "_key_gene\\.tsv$",
                          full.names = TRUE)
events <- data.frame(
  event_category = character(), gene = character(), celltype = character(),
  location = character(), patient = character(), cancer_type = character()
)
for (i in event_files) {
  x <- read.delim(i, stringsAsFactors = FALSE)
  if (nrow(x) == 0) next
  x$location <- ifelse(x$event_category == "StageSpecific", x$stage,
                       paste(x$from_stage, x$to_stage, sep = "->"))
  events <- rbind(events, x[, c("event_category", "gene", "celltype",
                                 "location", "patient", "cancer_type")])
}
events <- unique(events)

prop_scores <- data.frame(
  cancer_type = character(), gene = character(), celltype = character(),
  location = character(), event_category = character(),
  recurrence_count = integer(), total_patient_analyses = integer(),
  recurrence_frequency = numeric(), PropScore = numeric()
)
for (category in c("StageSpecific", "CrossStageConsistent")) {
  if (category == "StageSpecific") {
    eligible <- coverage[coverage$StageSpecific, ]
  } else {
    eligible <- coverage[coverage$CrossStageConsistent, ]
  }
  x <- events[events$event_category == category, ]
  keys <- unique(x[, c("cancer_type", "gene", "celltype", "location")])
  if (nrow(keys) == 0) next

  keys$event_category <- category
  keys$recurrence_count <- 0
  keys$total_patient_analyses <- 0
  for (j in seq_len(nrow(keys))) {
    same_event <- x$cancer_type == keys$cancer_type[j] &
      x$gene == keys$gene[j] & x$celltype == keys$celltype[j] &
      x$location == keys$location[j]
    same_cancer <- eligible$cancer_type == keys$cancer_type[j]
    keys$recurrence_count[j] <-
      length(intersect(unique(x$patient[same_event]), eligible$patient[same_cancer]))
    keys$total_patient_analyses[j] <- sum(same_cancer)
  }
  keys$recurrence_frequency <- keys$recurrence_count /
    keys$total_patient_analyses
  keys$PropScore <- NA_real_
  for (cancer in unique(keys$cancer_type)) {
    selected <- keys$cancer_type == cancer
    spread <- sd(keys$recurrence_frequency[selected])
    if (sum(selected) > 1 && is.finite(spread) && spread > 0) {
      keys$PropScore[selected] <- as.numeric(scale(keys$recurrence_frequency[selected]))
    }
  }
  prop_scores <- rbind(prop_scores, keys)
}

write.table(prop_scores,
            "outputs/database_events/molecular_event_PropScore.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

gene_files <- list.files("outputs/trajectory_gene",
                         pattern = "_trajectory_gene\\.tsv$",
                         full.names = TRUE)
trajectory_names <- character()
trajectory_genes <- list()

for (i in gene_files) {
  x <- read.delim(i, stringsAsFactors = FALSE)
  genes <- unique(x$gene[x$trajectory_associated %in% TRUE])
  if (length(genes) == 0) next
  trajectory_names <- c(trajectory_names,
                        sub("_trajectory_gene\\.tsv$", "", basename(i)))
  trajectory_genes[[length(trajectory_genes) + 1]] <- genes
}

jaccard <- data.frame(trajectory_1 = character(), trajectory_2 = character(),
                      shared_genes = integer(), jaccard = numeric())
if (length(trajectory_genes) > 1) {
  for (i in seq_len(length(trajectory_genes) - 1)) {
    for (j in (i + 1):length(trajectory_genes)) {
      shared <- length(intersect(trajectory_genes[[i]], trajectory_genes[[j]]))
      score <- shared / length(union(trajectory_genes[[i]], trajectory_genes[[j]]))
      if (score > 0.4 && score < 1) {
        jaccard <- rbind(jaccard, data.frame(
          trajectory_1 = trajectory_names[i],
          trajectory_2 = trajectory_names[j],
          shared_genes = shared,
          jaccard = score
        ))
      }
    }
  }
}
write.table(jaccard,
            "outputs/trajectory_similarity/trajectory_jaccard.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Step 5 complete.")
