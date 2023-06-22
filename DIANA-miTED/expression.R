# Read in the miRNA table
miRNA_table <- read.table("miRNA_table.csv", header = TRUE, sep = "\t", check.names = FALSE)

# Identify the columns containing miRNA expression data
miRNA_cols <- grep("^hsa-.", colnames(miRNA_table))

# Subset the miRNA expression table to include only the miRNA columns
miRNA_expr <- miRNA_table[, miRNA_cols]

# Calculate the median count value for each miRNA
median_counts <- apply(miRNA_expr, 2, median)

# Filter out miRNAs with low expression (count < 25th percentile value for that miRNA)
quantile_counts <- apply(miRNA_expr, 2, quantile, probs = 0.25)
miRNA_expr_filt <- miRNA_expr[, median_counts > quantile_counts[2]]

# output 
write.table(colnames(miRNA_expr_filt), file = 'miRNA_expr_filt.txt', row.names = F)

# Get the number of miRNAs retained after filtering
num_miRNAs_retained <- dim(miRNA_expr_filt)[2]
cat(sprintf("Number of miRNAs retained after filtering: %d\n", num_miRNAs_retained))
