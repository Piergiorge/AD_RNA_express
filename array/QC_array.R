library(pvca)
library(limma)
library(affy)
library(lme4)
library(PCAtools)
library(arrayQualityMetrics)
library(ggplot2)
library(gridExtra)

# Build Function to Return Element Text Object
rotatedAxisElementText <- function(angle, position = 'x') {
  angle <- angle[1]
  position <- position[1]
  positions <- list(x = 0, y = 90, top = 180, right = 270)
  
  if (!position %in% names(positions)) {
    stop(sprintf("'position' must be one of [%s]", paste(names(positions), collapse = ", ")), call. = FALSE)
  }
  
  if (!is.numeric(angle)) {
    stop("'angle' must be numeric", call. = FALSE)
  }
  
  rads <- (angle - positions[[position]]) * pi / 180
  hjust <- 0.5 * (1 - sin(rads))
  vjust <- 0.5 * (1 + cos(rads))
  element_text(angle = angle, vjust = vjust, hjust = hjust)
}

setwd("GSE5281/")
datadir <- paste(c(getwd(), "/"), collapse = "")

# Metadata
targets <- read.table("pdata.txt", sep = "\t", header = TRUE, row.names = 1)
targets[targets$Disease %in% c("normal", "Normal"), ]$Disease <- "Control"
targets[targets$Sex %in% c("male", "m", "M"), ]$Sex <- "Male"
targets[targets$Sex %in% c("female", "f", "F"), ]$Sex <- "Female"
targets[targets$Disease %in% c("Incipient", "Moderate", "Severe"), ]$Disease <- "AD"
targets[targets$Disease %in% c("definite AD", "Probable AD"), ]$Disease <- "AD"
targets[targets$Date == "11/01/04", ]$Date <- "2004"  # grouping with the nearest date
targets[targets$Date == "10/28/04", ]$Date <- "2004"
targets[targets$Date == "07/27/05", ]$Date <- "2005"
targets[targets$Date == "06/02/06", ]$Date <- "2006"  # only one on 06/02/06
targets[targets$Date == "07/10/06", ]$Date <- "2006"
targets

Data <- ReadAffy(filenames = row.names(targets), celfile.path = datadir, phenoData = targets)
pData(Data)

### PVCA cutoff
pct_threshold <- 0.1

# Adapt using the data documented in the metadata.
batch.factors <- c("Age", "Disease", "Date", "Sex", "Cell.Type", "Ethnicity")
pvcaObj <- pvcaBatchAssess(Data, batch.factors, pct_threshold)

dataf <- data.frame(pvcaObj$dat[, ], pvcaObj$label)
colnames(dataf) <- c('dat', 'label')
dataf

svg(filename = "raw_pvca.svg")
ggplot(data = dataf, aes(y = dataf$dat, x = reorder(dataf$label, -dataf$dat))) +
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab('Groups') +
  ylab("Weighted average proportion variance") +
  ggtitle("Raw data") +
  ylim(c(0, 1)) +
  theme(axis.text.x = rotatedAxisElementText(45, 'x'), axis.text.y = rotatedAxisElementText(45, 'y'))
dev.off()

preparedData <- prepdata(
  expressionset = Data,
  intgroup = c(),
  do.logtransform = TRUE
)

bo <- aqm.boxplot(preparedData)
de <- aqm.density(preparedData)
qpca <- aqm.pca(preparedData)
qheat <- aqm.heatmap(preparedData)
qmaplt <- aqm.maplot(preparedData)
qmean <- aqm.meansd(preparedData)
qm <- list("Boxplot" = bo, "Density" = de, "PCA" = qpca, "HeatMap" = qheat, "MAplot" = qmaplt,
           "Mean" = qmean)

bo@outliers@which
de@outliers@which
qpca@outliers@which
qheat@outliers@which
qmaplt@outliers@which
qmean@outliers@which

outdir <- paste(c(datadir, "arrayQualityMetrics_raw"), collapse = "")
aqm.writereport(
  modules = qm,
  reporttitle = "arrayQualityMetrics",
  outdir = outdir,
  arrayTable = pData(Data)
)
outdir


### normalization using RMA
eset <- rma(Data)
batch.factors <- c("Age", "Disease", "Date", "Sex", "Cell.Type", "Ethnicity")
pvcaObj <- pvcaBatchAssess (eset, batch.factors, pct_threshold)
pvcaObj

dataf <- data.frame(pvcaObj$dat[,], pvcaObj$label)
colnames(dataf) <- c('dat', 'label')
dataf

svg(filename="norm_pvca.svg")
ggplot(data = dataf, aes(y = dataf$dat, x = reorder(dataf$label, -dataf$dat))) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  xlab('Groups') + 
  ylab("Weighted average proportion variance") + 
  ggtitle("Normalized data") +
  ylim(c(0,1)) +
  theme(axis.text.x = rotatedAxisElementText(45,'x'), axis.text.y = rotatedAxisElementText(45,'y'))
dev.off()

### intergroup = 	the name of the sample covariate(s) used to draw a colour sidebar next to the heatmap.
### The first element of intgroup is also used to define sample groups in other plots (boxplots, densities).
### intgroup should be a character vector, and its elements must match the column's names of phenoData(expressionset).
### If its length is 0, the plots are not decorated with sample covariate information.
preparedData = prepdata(expressionset = eset,
                        intgroup = c(),
                        do.logtransform = TRUE)

bo = aqm.boxplot(preparedData)
de = aqm.density(preparedData)
qpca = aqm.pca (preparedData)
qheat = aqm.heatmap (preparedData)
qmaplt = aqm.maplot(preparedData)
qmean = aqm.meansd(preparedData)

qm = list("Boxplot" = bo, "Density" = de, "PCA" = qpca, "HeatMap" = qheat, "MAplot" =  qmaplt, 
          "Mean" = qmean)
#bo@outliers
bo@outliers@which
de@outliers@which
qpca@outliers@which
qheat@outliers@which
qmaplt@outliers@which
qmean@outliers@which

outdir = paste (c(datadir,"arrayQualityMetrics_normalization"), collapse = "")
aqm.writereport(modules = c(qm), reporttitle = "arrayQualityMetrics", outdir = outdir,
                arrayTable = pData(eset))
outdir

### PCA - raw 
ex <- exprs(Data)
# Remove this % of variables based on low variance.
princip <- pca(ex, metadata = targets, removeVar = 0.1)

svg(filename = "plot_bi1.svg", width = 14, height = 7)
  biplot(princip, showLoadings = TRUE, labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
dev.off()

svg(filename = "plot1.svg", width = 14, height = 7)
  screeplot(princip, axisLabSize = 18, titleLabSize = 22)
dev.off()

svg(filename = "plot2.svg", width = 14, height = 7)
  biplot(princip, lab = Data$Age, hline = 0, vline = 0, legendPosition = 'right')
dev.off()

svg(filename = "plot3.svg", width = 14, height = 7)
  biplot(princip, lab = Data$Sex, hline = 0, vline = 0, legendPosition = 'right')
dev.off()

svg(filename = "plot4.svg", width = 14, height = 7)
  biplot(princip, lab = Data$Disease, hline = 0, vline = 0, legendPosition = 'right')
dev.off()

svg(filename = "plot5.svg", width = 14, height = 7)
  biplot(princip, lab = Data$Date, hline = 0, vline = 0, legendPosition = 'right')
dev.off()

### PCA - RMA 
ex <- exprs(eset)
# Remove this % of variables based on low variance.
princip <- pca(ex, metadata = targets, removeVar = 0.1)

svg(filename = "plot_bi1_rma.svg", width = 14, height = 7)
  biplot(princip, showLoadings = TRUE, labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
dev.off()

svg(filename = "plot1_rma.svg", width = 14, height = 7)
  screeplot(princip, axisLabSize = 18, titleLabSize = 22)
dev.off()

svg(filename = "plot2_rma.svg", width = 14, height = 7)
  biplot(princip, lab = eset$Age, hline = 0, vline = 0, legendPosition = 'right')
dev.off()

svg(filename = "plot3_rma.svg", width = 14, height = 7)
  biplot(princip, lab = eset$Sex, hline = 0, vline = 0, legendPosition = 'right')
dev.off()

svg(filename = "plot4_rma.svg", width = 14, height = 7)
  biplot(princip, lab = eset$Disease, hline = 0, vline = 0, legendPosition = 'right')
dev.off()

svg(filename = "plot5_rma.svg", width = 14, height = 7)
  biplot(princip, lab = eset$Date, hline = 0, vline = 0, legendPosition = 'right')
dev.off()

### Plot the density distributions before and after normalization
svg(filename = "density_before.svg", width = 14, height = 7)
  plotDensities(log2(exprs(Data)), main = "Before normalization (log2-transformed) - Probe level", legend = "topright")
dev.off()

svg(filename = "density_after.svg", width = 14, height = 7)
  plotDensities(exprs(eset), main = "RMA-normalized - Probeset level", legend = "topright")
dev.off()
