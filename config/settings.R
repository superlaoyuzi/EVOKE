# Parameters stated in the manuscript Methods or used for the Pt-10 example.

STAGE_ORDER <- c("Health", "Inflamed", "Precancerous", "Margin",
                 "Primary", "Metastasis")

MIN_FEATURES <- 200
MAX_MITO_PERCENT <- 20
NORMALIZATION_SCALE <- 10000
N_VARIABLE_FEATURES <- 2000
N_PCS <- 15
CLUSTER_RESOLUTION <- 2.0
RANDOM_SEED <- 10086

SCTYPE_WRAPPER <- "sctype/sctype_wrapper.R"
SCTYPE_MARKERS <- "sctype/ScTypeDB_full.xlsx"
SCTYPE_TISSUE <- "Stomach"

AUGUR_HIGH_CONFIDENCE_CUTOFF <- 0.8
ROE_ENRICHMENT_CUTOFF <- 2
ROE_DEPLETION_CUTOFF <- 0.5

TRADESEQ_MIN_CELL_FRACTION <- 0.05
TRADESEQ_NKNOTS <- 6
TRADESEQ_WORKERS <- 2L # Set to 1L if memory is limited.
TRADESEQ_RUN_KNOT_CHECK <- TRUE # Set FALSE to skip the optional evaluateK diagnostic.
TRADESEQ_P_CUTOFF <- 0.05
TRADESEQ_LOG2FC_CUTOFF <- 1
GO_BP_ADJ_P_CUTOFF <- 0.01
