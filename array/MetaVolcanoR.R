library(MetaVolcanoR)

# Alzheimer
GSE5281 <- read.table('results/total_res_GSE5281.txt', sep = '\t', header = TRUE)
GSE48350 <- read.table('results/total_res_GSE48350.txt', sep = '\t', header = TRUE)
GSE28146 <- read.table('results/total_res_GSE28146.txt', sep = '\t', header = TRUE)
GSE1297 <- read.table('results/total_res_GSE1297.txt', sep = '\t', header = TRUE)
GSE36980 <- read.table('results/total_res_GSE36980.txt', sep = '\t', header = TRUE)
GSE84422 <- read.table('results/total_res_GSE84422.txt', sep = '\t', header = TRUE)

dataset <- list()
dataset[['GSE5281']] <- GSE5281
dataset[['GSE48350']] <- GSE48350
dataset[['GSE28146']] <- GSE28146
dataset[['GSE1297']] <- GSE1297
dataset[['GSE36980']] <- GSE36980
dataset[['GSE84422']] <- GSE84422

# The REM MetaVolcano summarizes the gene fold change of several studies taking into account the variance.
# The REM estimates a summary p-value which stands for the probability of the summary fold-change is not different than zero.
# Users can set the metathr parameter to highlight the top percentage of the most consistently perturbed genes.
# This perturbation ranking is defined following the topconfects approach.

meta_degs_rem <- rem_mv(
  diffexp = dataset,
  pcriteria = "pvalue", #the column name of the pvalue variable <string>
  foldchangecol = 'log2FC',  #the column name of the foldchange variable <string>
  genenamecol = 'Symbol', #the column name of the genename variable <string>
  geneidcol = NULL,
  collaps = FALSE,
  llcol = 'CI.L',
  rlcol = 'CI.R',
  vcol = NULL,
  cvar = TRUE,
  metathr = 0.01, #top percentage of perturbed genes to be highlighted <double> 0.01
  jobname = "MetaVolcano",
  outputfolder = ".",
  draw = 'PDF'
)

# REM results
head(meta_degs_rem@metaresult, 3)

write.table(meta_degs_rem@metaresult, file = "meta_degs_rem.tsv", row.names = FALSE, sep = "\t")

# Plot MetaVolcano
svg(filename = "meta_degs_rem_volcano.svg")
meta_degs_rem@MetaVolcano
dev.off()
