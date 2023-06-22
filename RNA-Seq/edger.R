library(edgeR)
library(ggplot2)

experiment <- "dir/pheno.txt"

# Getting experiment data.
design.table <- read.table(experiment, sep = "\t", header = TRUE)
archs <- as.vector(design.table[, 1])  # Column names
cond1 <- as.vector(design.table[, 2])  # Cell types

counts.table <- read.table("dir/counts.txt", header = TRUE, sep = '\t')
rownames(counts.table) <- paste(counts.table$ENSEMBL, counts.table$GENE, sep = '_')
counts.table$GENE <- NULL
counts.table$ENSEMBL <- NULL

counts.table <- counts.table

# Getting group condition;
group <- vector()
count.order <- noquote(as.vector(colnames(counts.table)))

for (d in 1:length(count.order)) {
  for (e in 1:length(archs)) {
    if (count.order[d] == archs[e]) {
      group[d] <- cond1[e]
    }
  }
}

group <- as.factor(group)

# Start differential expression analysis
rm(y)
y <- DGEList(counts = counts.table, group = group)
dim.table <- dim(y)
dim.table

keep <- rowSums(cpm(y) > 1) >= dim.table[2] / 2
y <- y[keep,]
dim(y)

# Recalculate library size
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
y <- estimateDisp(y)

et <- exactTest(y, pair = c("AD", "Control"))
de.top <- topTags(et, n = length(et$table[, 1]))
de.cpm <- cpm(y)[rownames(de.top),]

# Writing outputs
de.top$table$Gene <- rownames(de.top)
de.cpm <- as.data.frame(de.cpm)
de.cpm$Gene <- rownames(de.cpm)
de.all <- merge(de.cpm, de.top, by = "Gene", all = TRUE)
de.all2 <- de.all[order(de.all$PValue, de.all$FDR < 0.05) & abs(de.all$logFC) > 1, ]
results <- de.top$table

de.top$table <- de.top$table[abs(de.top$table$logFC) > 1 & de.top$table$FDR < 0.05, ]
write.table(de.top$table, file = "all_expression_control_vs_AD.txt", row.names = FALSE, sep = "\t")
