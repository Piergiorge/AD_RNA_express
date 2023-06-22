  # AD_RNA_express - Array

This folder contains the array-related files for the AD_RNA_express project.

## Description

The "array" folder within the AD_RNA_express project repository contains files and scripts related to array data processing and analysis. It includes raw data files, pre-processing scripts, quality control analyses, and visualisations.


# series_matrix2pdata.sh

This code is utilized to convert a GSE..._series_matrix.txt file into a TSV format. Although additional manual editing may be necessary, this code streamlines the sample handling process.

## Usage

```bash
./series_matrix2pdata.sh GSE20333_series_matrix.txt
```

# QC_array.R:
The "QC_array.R" file is an R script that implements quality control (QC) for microarray data. It utilizes R libraries and functions to perform the microarray data quality assessment. The script includes steps for data preprocessing, generating quality plots and metrics, detecting outliers, performing differential expression analysis, and visualizing the results.

# expression_affy.R:
The "expression_affy.R" file is an R script that analyses gene expression using Affymetrix microarray data. It utilizes the Affy library and other related R libraries to load the microarray data, preprocess it, perform normalization, identify differentially expressed genes, and visualize the results. The script includes steps like reading the microarray data files, gene filtering, data normalization, and conducting statistical tests to identify significantly expressed genes. It is a valuable tool for gene expression analysis from Affymetrix microarray data related to Alzheimer's disease.
