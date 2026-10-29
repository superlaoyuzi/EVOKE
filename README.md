# EVOKE: A Single-Cell Resource for Mapping Tumor Progression

EVOKE (Evolutionary Vector of Key Events) is a publicly accessible database for exploring progression-associated cellular and molecular changes in human and mouse cancers at single-cell resolution. This repository hosts the core analysis scripts associated with the EVOKE study.

## Overview

Tumor progression involves changes in cellular composition, transcriptional states, differentiation trajectories, and the tumor microenvironment. EVOKE organizes publicly available single-cell RNA-sequencing (scRNA-seq) datasets by annotated disease states and progression processes to facilitate systematic exploration of these changes.

The current EVOKE release contains:

| Resource | Current release |
| --- | ---: |
| Single-cell transcriptomes | **7,147,991** |
| Human and mouse samples | **1,141** |
| Cancer types | **12** |
| Defined progression axes | **8** |
| Progression-associated molecular events | **>29,000** |

Datasets were curated from public repositories, including [GEO](https://www.ncbi.nlm.nih.gov/geo/) and [ArrayExpress](https://www.ebi.ac.uk/biostudies/arrayexpress), with a search cutoff of **June 2026**. Inclusion was based on progression-relevant annotations, available expression matrices, sample metadata, and quality-control requirements. Dataset accessions and characteristics are listed in Supplementary Table S1 of the accompanying manuscript.

## Key Features

- **Progression-oriented data organization:** Explore annotated healthy, inflammatory, primary-tumor, and metastatic states across defined progression contexts.
- **Cell Events:** Characterize cell-type composition, transcriptional perturbation, and tissue preference during progression.
- **Trajectory Events:** Explore inferred cellular trajectories, pseudotime-dependent expression patterns, and lineage-associated changes.
- **Gene Events:** Examine stage-specific perturbation genes and cross-stage consistently altered genes, together with related functional annotations.
- **Interactive exploration:** Browse and search by cancer type, progression context, cell type, or gene; access analysis tools and downloadable results through the EVOKE website.

EVOKE characterizes *transcriptome-associated* progression patterns; inferred trajectories and molecular events should not be interpreted as direct evidence of genomic clonal evolution or experimentally established causal mechanisms.

## Repository Contents

The repository is organized around the following analysis steps and includes a three-stage example dataset:

| File or directory | Purpose |
| --- | --- |
| [`example_data.zip`](example_data.zip) | Processed Pt-10 gene-by-cell UMI count matrices for three disease states |
| [`config/settings.R`](config/settings.R) | Parameters used by the example workflow |
| [`scripts/00_prepare_pt10.R`](scripts/00_prepare_pt10.R) | Step 0: prepare example expression matrices and sample metadata |
| [`scripts/01_preprocess_and_annotate.R`](scripts/01_preprocess_and_annotate.R) | Step 1: quality control, preprocessing, and cell-type annotation |
| [`scripts/02_cellular_events.R`](scripts/02_cellular_events.R) | Step 2: cellular composition, tissue preference, and perturbation analysis |
| [`scripts/03_trajectory_events.R`](scripts/03_trajectory_events.R) | Step 3: Slingshot trajectory inference and pseudotime analysis |
| [`scripts/04_1_key_genes_in_progression.R`](scripts/04_1_key_genes_in_progression.R) | Step 4.1: progression-associated gene events |
| [`scripts/04_2_trajectory_associated_genes.R`](scripts/04_2_trajectory_associated_genes.R) | Step 4.2: trajectory-associated gene analysis |
| [`scripts/05_recurrence_and_trajectory_similarity.R`](scripts/05_recurrence_and_trajectory_similarity.R) | Step 5: recurrence scoring and trajectory similarity |
| [`scripts/06_GO_BP_and_AvgRichFactor.R`](scripts/06_GO_BP_and_AvgRichFactor.R) | Step 6: GO biological-process enrichment and summary metrics |
| [`scripts/run_all.R`](scripts/run_all.R) | Execute the example workflow in sequence |
| [`sctype/`](sctype/) | Supporting code and marker resources for scType annotation |

## Example Data and Usage

A small example from **GSE239676** is provided to demonstrate EVOKE's workflow using three matched tissue samples from a single patient (**Pt-10**) with gastric cancer (STAD):

| Sample ID | GEO accession | Tissue | EVOKE progression state | Cells |
| --- | --- | --- | --- | ---: |
| `Pt10_Ad` | `GSM7669805` | Adjacent normal gastric mucosa | Health | 1,200 |
| `Pt10_P` | `GSM7669807` | Primary gastric tumor | Primary | 1,200 |
| `Pt10_Li` | `GSM7669806` | Liver metastasis | Metastasis | 1,200 |

The supplied [`example_data.zip`](example_data.zip) contains **3,600 preselected cells** in three 10x-style gene-by-cell UMI-count matrices (`matrix.mtx.gz`, `features.tsv.gz`, and `barcodes.tsv.gz` for each sample). These are processed expression counts, **not raw sequencing FASTQ files**. The three samples represent different disease states of one patient, not three longitudinal collection time points. The subset is intended to demonstrate the workflow rather than reproduce the full EVOKE database results.

### Run the Pt-10 example

Clone or download this repository, open R, and set the working directory to the repository root (the directory containing `example_data.zip`). With the required R packages installed, run the complete example:

```r
source("scripts/run_all.R")
```

To run or inspect the individual steps instead:

```r
source("scripts/00_prepare_pt10.R")
source("scripts/01_preprocess_and_annotate.R")
source("scripts/02_cellular_events.R")
source("scripts/03_trajectory_events.R")
source("scripts/04_1_key_genes_in_progression.R")
source("scripts/04_2_trajectory_associated_genes.R")
source("scripts/05_recurrence_and_trajectory_similarity.R")
source("scripts/06_GO_BP_and_AvgRichFactor.R")
```

**Step 0** automatically reads and extracts `example_data.zip`, combines the three samples with their progression-state metadata, and saves `data/Pt10_merged_raw.rds` for **Step 1**. **Step 1** performs quality control and cell-type annotation and writes `outputs/stage1/STAD_Pt-10_after_cell_annotation.rds`. The subsequent steps save their results under `outputs/`. There is no need to unzip the example archive or enter sample paths manually.

For optional manual review of cell-type assignments, Step 1 also generates `outputs/stage1/STAD_Pt-10_cluster_celltype_review.csv`. The `reviewed_label` column can be filled only where a correction is supported by marker evidence; rerun Step 1 after making any changes. Blank labels retain the automatic annotation.

Example parameters are specified in [`config/settings.R`](config/settings.R). The workflow requires Seurat and its associated R dependencies, as well as packages used in later steps (including Augur, Slingshot, tradeSeq, and clusterProfiler) and the bundled scType resources. The example is not a substitute for the complete multi-patient, multi-cancer production analysis: in particular, a single patient is insufficient to estimate meaningful cross-patient recurrence scores. The packaged workflow has not yet been independently validated end to end in an R environment; please verify the generated outputs before using it for reproduction claims.

## Data Processing and Reproducibility

EVOKE applies a standardized, sample-wise analysis strategy. Samples are processed individually, and the resulting annotations and progression-associated events are harmonized at the database level by cancer type, species, progression state, cell type, dataset, and sample identity. This database-level harmonization should not be confused with merging cells from unrelated studies into one corrected expression matrix.

For the bundled Pt-10 example, the expression matrices and matched sample metadata are already supplied; Step 0 prepares the input for downstream analysis. For other datasets, prepare gene-by-cell expression matrices and matched metadata, including sample identifiers, species, tissue source, cancer type, and progression-state annotations, and update the relevant configuration and input paths. The example workflow runs from data preparation (Step 0) to preprocessing (Step 1), cellular events (Step 2), trajectory analysis (Step 3), and molecular analyses (Steps 4–6).

For source-data accessions, preprocessing criteria, statistical definitions, software versions, and interpretation of outputs, refer to the **Methods** and **Supplementary Materials** of the associated EVOKE manuscript. The analysis scripts should be used with the finalized production settings corresponding to the reported release.

## Web Access

The interactive EVOKE database is available at:

**[https://bio-bigdata.hrbmu.edu.cn/evoke/](https://bio-bigdata.hrbmu.edu.cn/evoke/)**

The website provides **Browse**, **Search**, **Tools**, and **Download** functions for exploration of progression-associated single-cell data and analysis results.

## Citation

If you use EVOKE or its analysis code, please cite the associated manuscript:

> Li X, Bai J, Li R, et al. **EVOKE: a comprehensive single-cell visualization resource for mapping evolutionary events in tumor progression.** *Genomics, Proteomics & Bioinformatics* (manuscript under revision).

Publication details and a DOI will be added when available.
