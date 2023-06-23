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

# MetaVolcanoR.R
MetaVolcanoR is an R script for performing a meta-analysis of gene expression data using the REM (Rank-based Empirical Bayes) MetaVolcano approach. It allows users to summarize gene fold change across multiple studies, estimate summary p-values, and identify consistently perturbed genes.

* <div class="csl-entry">Prada-Medina, C. A., Peron, J. P. S., &#38; Nakaya, H. I. (2020). Immature neutrophil signature associated with the sexual dimorphism of systemic juvenile idiopathic arthritis. <i>Journal of Leukocyte Biology</i>, <i>108</i>(4), 1319–1327. https://doi.org/10.1002/JLB.6MA0720-015RR</div>

# integrate.R

The goal of the integrate.R script is to merge RNA expression data from various microarray datasets associated with AD. The process involves reading and processing the data, utilizing the CEMiTool package for module detection, computing overlaps between modules, extracting significant gene information, and generating a file called `cemoverlap_df.txt`. This file contains information on module overlaps, including gene pairs and edge counts.


* <div class="csl-entry">Russo, P. S. T., Ferreira, G. R., Cardozo, L. E., Bürger, M. C., Arias-Carrasco, R., Maruyama, S. R., Hirata, T. D. C., Lima, D. S., Passos, F. M., Fukutani, K. F., Lever, M., Silva, J. S., Maracaja-Coutinho, V., &#38; Nakaya, H. I. (2018). CEMiTool: a Bioconductor package for performing comprehensive modular co-expression analyses. <i>BMC Bioinformatics</i>, <i>19</i>(1), 56. https://doi.org/10.1186/s12859-018-2053-1</div>

* <div class="csl-entry">de Lima, D. S., Cardozo, L. E., Maracaja-Coutinho, V., Suhrbier, A., Mane, K., Jeffries, D., Silveira, E. L. v, Amaral, P. P., Rappuoli, R., de Silva, T. I., &#38; Nakaya, H. I. (2019). Long noncoding RNAs are involved in multiple immunological pathways in response to vaccination. <i>Proceedings of the National Academy of Sciences of the United States of America</i>, <i>116</i>(34), 17121–17126. https://doi.org/10.1073/pnas.1822046116</div>
