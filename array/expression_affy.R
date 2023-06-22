library(affy)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(limma)
library(lme4)
library(stringr)

### Set Path
setwd("GSE48350/")
datadir <- paste(c(getwd(), "/"), collapse = "")

### Metadata
targets <- read.table("pdata.txt", sep = "\t", header = TRUE, row.names = 1)
targets[targets$Disease == "normal" | targets$Disease == "Normal", ]$Disease <- "Control"
targets[targets$Sex == "male" | targets$Sex == "m" | targets$Sex == "M", ]$Sex <- "Male"
targets[targets$Sex == "female" | targets$Sex == "f" | targets$Sex == "F", ]$Sex <- "Female"
targets[targets$Disease %in% c("Incipient", "Moderate", "Severe"), ]$Disease <- "AD"
targets[targets$Disease %in% c("definite AD", "Probable AD"), ]$Disease <- "AD"
targets[targets$Date == "11/30/05", ]$Date <- "2005" # grouping with nearby date
targets[targets$Date == "12/02/05", ]$Date <- "2005"
targets[targets$Date == "12/05/05", ]$Date <- "2005"
targets[targets$Date == "10/23/07", ]$Date <- "2007"
targets[targets$Date == "10/24/07", ]$Date <- "2007"
targets[targets$Date == "02/08/06", ]$Date <- "2006a"
targets[targets$Date == "02/02/06", ]$Date <- "2006a"
targets[targets$Date == "02/03/06", ]$Date <- "2006a"
targets[targets$Date == "01/18/06", ]$Date <- "2006b"
targets[targets$Date == "01/13/06", ]$Date <- "2006b"
targets[targets$Date == "01/10/06", ]$Date <- "2006b"
targets[targets$Date == "06/22/06", ]$Date <- "2006c"
targets[targets$Date == "06/23/06", ]$Date <- "2006c"

Data <- ReadAffy(filenames = row.names(targets), celfile.path = datadir, phenoData = targets)

### Normalization (RMA)
eset <- rma(Data)
ex <- exprs(eset)

# Filter > 25% & > min sample size
quantil25 <- quantile(ex)[2]
no_of_samples <- table(eset$Disease)
samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(exprs(eset), 1, function(x) sum(x > quantil25) >= samples_cutoff)
table(idx_man_threshold)
Data_manfiltered_1 <- subset(exprs(eset), idx_man_threshold)

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annot_conf <- getBM(
  mart = mart,
  attributes = c(
    "affy_hg_u133_plus_2",
    "chromosome_name",
    "ensembl_gene_id",
    "external_gene_name"
  ),
  filter = "affy_hg_u133_plus_2",
  values = rownames(Data_manfiltered_1),
  uniqueRows = TRUE
)
dim(annot_conf)

annot_conf <- as.data.frame(annot_conf)
names(annot_conf)[1] <- "name"

# Assign Ensembl ID to probes without gene name
annot_conf[annot_conf == ''] <- NA
annot_conf <- annot_conf %>% mutate(across(everything(), na_if, '')) %>%
  mutate(external_gene_name = coalesce(!!!select(., external_gene_name:ensembl_gene_id)))
chrom <- c(1:22, 'X', 'Y', 'MT')
annot_conf <- annot_conf[annot_conf$chromosome_name %in% chrom, ]

temp <- as.data.frame(Data_manfiltered_1)
temp$name <- rownames(temp)
temp <- merge(temp, annot_conf, by = 'name')

data_summarised <- as.data.frame(avereps(temp, ID = temp$external_gene_name))
rownames(data_summarised) <- data_summarised$external_gene_name
data_summarised$ensembl_gene_id <- NULL
data_summarised$chromosome_name <- NULL
data_summarised$external_gene_name <- NULL
data_summarised$name <- NULL

# Create matrix for analysis
colunas <- colnames(data_summarised)
linhas <- rownames(data_summarised)
matriz <- sapply(data_summarised, as.numeric)
dimnames(matriz) <- list(linhas, colunas)

arquivo <- str_match(datadir, "GSE[0-9]+")
arquivo <- paste(c("matriz", "_", arquivo, ".txt"), collapse = "")
write.table(matriz, file = arquivo, sep = "\t", row.names = TRUE)

### Metadata for analysis
cond <- factor(eset$Disease)
cond <- relevel(x = cond, ref = "AD")
Date <- factor(eset$Date)
Sex <- factor(eset$Sex)
Age <- factor(eset$Age)
design <- model.matrix(~0 + cond + as.numeric(Sex) + as.numeric(Date))
colnames(design) <- c('AD', 'Control', 'Sex', 'Date')
fit <- lmFit(matriz, design) # Fit linear model
con <- makeContrasts("AD-Control", levels = design)
fit2 <- contrasts.fit(fit, con)
fit2 <- eBayes(fit2)
options(digits = 2)
total_res <- topTreat(fit2, coef = 1, number = Inf, adjust.method = "fdr", confint = TRUE, sort.by = "p")
head(total_res)

# Output
diffexp <- data.frame(Symbol = rownames(total_res),
                      log2FC = total_res$logFC,
                      pvalue = total_res$P.Value,
                      CI.L = total_res$CI.L,
                      CI.R = total_res$CI.R)

arquivo <- str_match(datadir, "GSE[0-9]+")
arquivo <- paste(c("total_res", "_", arquivo, ".txt"), collapse = "")
write.table(diffexp, file = arquivo, row.names = FALSE, sep = "\t")
