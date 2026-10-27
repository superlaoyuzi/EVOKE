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

The repository is organized around the principal analysis steps:

| File or directory | Purpose |
| --- | --- |
| [`Step1_Data preprocessing and celltype annotation.R`](Step1_Data%20preprocessing%20and%20celltype%20annotation.R) | Sample-level preprocessing and cell-type annotation |
| [`Step2_Cellular events.R`](Step2_Cellular%20events.R) | Identification of progression-associated cellular events |
| [`Step3_Trajectory events.R`](Step3_Trajectory%20events.R) | Trajectory reconstruction and pseudotime analysis |
| [`Step4_1_keygene_in_progression.R`](Step4_1_keygene_in_progression.R) | Identification of progression-associated gene events across stages |
| [`Step4_2_trajectory_gene.R`](Step4_2_trajectory_gene.R) | Analysis of gene-expression dynamics along inferred trajectories |
| [`sctype/`](sctype/) | Supporting files for cell-type annotation |

## Data Processing and Reproducibility

EVOKE applies a standardized, sample-wise analysis strategy. Samples are processed individually, and the resulting annotations and progression-associated events are harmonized at the database level by cancer type, species, progression state, cell type, dataset, and sample identity. This database-level harmonization should not be confused with merging cells from unrelated studies into one corrected expression matrix.

To use the analysis scripts, prepare the relevant gene-by-cell expression matrices and matched sample-level metadata, including sample identifiers, species, tissue source, cancer type, and available progression-state annotations. Configure dataset paths and other required inputs in the corresponding scripts before execution. The workflow is organized from preprocessing and annotation (Step 1) to cellular events (Step 2), trajectory analysis (Step 3), and gene-event analyses (Steps 4.1 and 4.2).

For source-data accessions, preprocessing criteria, statistical definitions, software versions, and interpretation of outputs, refer to the **Methods** and **Supplementary Materials** of the associated EVOKE manuscript. The analysis scripts should be used with the finalized production settings corresponding to the reported release.

## Web Access

The interactive EVOKE database is available at:

**[https://bio-bigdata.hrbmu.edu.cn/evoke/](https://bio-bigdata.hrbmu.edu.cn/evoke/)**

The website provides **Browse**, **Search**, **Tools**, and **Download** functions for exploration of progression-associated single-cell data and analysis results.

## Citation

If you use EVOKE or its analysis code, please cite the associated manuscript:

> Li X, Bai J, Li R, et al. **EVOKE: a comprehensive single-cell visualization resource for mapping evolutionary events in tumor progression.** *Genomics, Proteomics & Bioinformatics* (manuscript under revision).

Publication details and a DOI will be added when available.
