# R-Project
PBMC Cell Analysis in R: From Gene Data to Machine Learning
# 🧬 Single-Cell RNA-Seq Analysis of PBMC 1k Dataset (Seurat + ML Classification)

**Objective:**  
Analyze single-cell RNA sequencing data to identify cell types and differentially expressed genes.
This repository contains an R script and dataset for analyzing single-cell RNA sequencing (scRNA-seq) data using **Seurat**, followed by **machine learning-based cell classification** using **caret** and **kNN**.

## 📂 Repository

**GitHub**: https://github.com/toobamir/scRNAseq-PBMC-Seurat-ML

## 📁 Contents

- `pbmc_analysis.R` — Main R script for the full Seurat workflow + machine learning
- `pbmc_1k_v2_filtered_feature_bc_matrix.h5` — Filtered feature-barcode matrix in HDF5 format from 10x Genomics
- `README.md` — Documentation and usage instructions

> **Note**: If the `.h5` file is too large for GitHub, it may not be uploaded. Download it from [10x Genomics PBMC 1k v2](https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_filtered_feature_bc_matrix.h5) and place it in your working directory.

---

## 🔧 Requirements

Install the following R packages before running the script:

```r
install.packages(c("Seurat", "patchwork", "dplyr", "ggplot2", "caret", "e1071", "randomForest"))
