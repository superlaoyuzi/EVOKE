# EVOKE: Single-cell Atlas of Cancer Evolution

This repository contains analysis code and processing scripts associated with the EVOKE resource, which systematically characterizes tumor progression using large-scale single-cell transcriptomic data.

## Overview

Cancer progression involves continuous evolutionary changes, spanning from normal or chronically inflamed tissues to primary tumors and eventually metastatic disease. These transitions are accompanied by extensive alterations in cellular composition, gene regulation, and microenvironmental interactions.

EVOKE is an integrated single-cell framework designed to systematically describe these dynamic processes across cancer evolution. It aggregates over 5.2 million single-cell transcriptomes derived from 797 patient samples, covering 12 cancer types across both human and mouse datasets. All datasets included are annotated with tumor progression stages, enabling a unified view of cancer development across conditions.

## Key Features

- **Large-scale integration of single-cell data** across multiple cancer types and species  
- **Stage-resolved analysis** from normal/inflammatory states to primary tumors and metastatic lesions  
- **Cellular state dynamics** including differentiation trajectories and compositional shifts  
- **Functional programs** such as gene expression signatures and pathway activity changes  
- **Evolutionary axes analysis** capturing multiple dimensions of tumor progression  
- Identification of thousands of **molecular events associated with malignant transformation and metastasis**

## Resource Description

The EVOKE framework reconstructs a continuous cellular and regulatory landscape of cancer evolution, enabling the exploration of:

- Microenvironmental remodeling during tumor progression  
- Lineage-specific differentiation dynamics  
- Regulatory network alterations driving malignancy  
- System-level changes associated with metastasis

In total, the dataset supports the characterization of tens of thousands of molecular events linked to cancer evolution.

## Repository Contents

This repository includes scripts for:

- Data preprocessing and integration  
- Downstream single-cell analysis  
- Trajectory and state inference workflows  
- Visualization of evolutionary patterns and functional programs  

## Web Access

An interactive version of the EVOKE resource is available at:  
:contentReference[oaicite:0]{index=0}

## Notes

This codebase is intended to support reproducibility of the analyses described in the associated study and to facilitate further exploration of cancer evolution using single-cell data.

## Citation

If you use this resource or code, please cite the original publication associated with EVOKE.
