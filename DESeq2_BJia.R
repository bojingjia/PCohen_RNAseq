# Initiate Bioconductor
source("http://bioconductor.org/biocLite.R")

# Load Libraries and Import Data
library('DESeq2')
library("gplots")
library("RColorBrewer")

# Import data from featureCounts 
countdata <- read.table("C:/Users/user/Downloads/Cohen Lab/RNAseq/E_GEOD_64459/IngVsEpi_counts.matrix", header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]

# Reorder the columns
# countdata <- countdata[,c(1,2,4,7,8,3,5,6,9,10)]
# Note[Feb 05 2016]: PCA shows batch effect, keep Ctrl 162, 161, 163 and KO 128, 139, 129
# countdata <- countdata[,c(1,4,7,6,9,10)]

# Convert to matrix
head(countdata)
countdata <- as.matrix(countdata)

# Assign sample condition to dataframe (an attribute that describes each assay )
#change name of "sampleCondition" to specify design formula
#sampleCondition<-c("Ing Ctrl", "Ing Ctrl","Ing Ctrl","Ing Ctrl","Ing Ctrl","Ing KO","Ing KO","Ing KO","Ing KO","Ing KO")
sampleDepot<-c("EWAT", "IWAT", "EWAT", "IWAT", "EWAT", "EWAT", "IWAT", "EWAT", "IWAT", "EWAT", "IWAT", "EWAT", "EWAT", "IWAT", "IWAT", "EWAT", "IWAT")
sampleGenotype<-c("129S", "129S", "B6", "B6", "129S", "129S", "129S", "B6", "B6", "B6", "129S", "129S", "B6", "B6", "B6", "B6", "129S")

(coldata <- data.frame(row.names=colnames(countdata), sampleDepot, sampleGenotype))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~sampleDepot + sampleGenotype + sampleDepot:sampleGenotype)

# To specify which condition is our reference level, we relevel our data
dds$sampleGenotype <- relevel(dds$sampleGenotype, "B6")
dds$sampleDepot <- relevel (dds$sampleDepot, "EWAT")

# Normalize read counts
dds<-DESeq(dds) #was previously DESeq(ddsHtseq), although it is unclear if Htseq is an option for the command

# First save normalized readcounts
setwd("C:/Users/user/R_sandbox")
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file = "IngVsEpi_B6v129S1_20161115_normalizedcounts.csv")

# Print out interaction terms
resultsNames(dds)

#Make binary contrasts for 2-factor experiment
res_IvsE_B6 <-results(dds, contrast=c("sampleDepot", "IWAT", "EWAT"))
res_IvsE_129S <-results(dds, contrast=list(c("sampleDepot_IWAT_vs_EWAT","sampleDepotIWAT.sampleGenotype129S")))
res_129SvsB6_E <-results(dds, contrast=c("sampleGenotype", "129S", "B6"))
res_129SvsB6_I <-results(dds, contrast=list(c("sampleGenotype_129S_vs_B6","sampleDepotIWAT.sampleGenotype129S")))
### Inherently, the interaction term is answering: is the condition effect *different* across genotypes?
### ie results(dds, name="genotypeII.conditionB")

# res<-res[order(res$padj),] to order by adjusted p-value
res_IvsE_B6 <-res_IvsE_B6[order(res_IvsE_B6$padj),]
res_IvsE_129S <- res_IvsE_129S[order(res_IvsE_129S$padj),]
res_129SvsB6_E <- res_129SvsB6_E[order(res_129SvsB6_E$padj),]
res_129SvsB6_I <- res_129SvsB6_I[order(res_129SvsB6_I$padj),]

# Write to spreadsheet
write.csv(res_IvsE_B6, file = "res_IvsE_B6.csv")
write.csv(res_IvsE_129S, file = "res_IvsE_129S.csv")
write.csv(res_129SvsB6_E, file = "res_129SvsB6_E.csv")
write.csv(res_129SvsB6_I, file = "res_129SvsB6_I.csv")

# for PCA, adds labels
##note: blind = FALSE means diff between replicates and treatment should not add to variance profile
rld <- rlogTransformation(dds, blind=FALSE)
d <- plotPCA(rld, intgroup =c("sampleDepot"), returnData=TRUE)
library(ggplot2)
ggplot(d, aes(x=PC1,y=PC2,col=sampleCondition,label=name)) + geom_point() + geom_text(nudge_x = 15, nudge_y = -0.2, check_overlap=TRUE, size = 3) + scale_x_continuous(expand = c(.1, 12.5))
dev.copy(png,"deseq2_pca_improved_labels.png")
dev.off()

# Find upregulated, and of those, significant genes
posLFC<-res[which(results_condition$log2FoldChange > 0),]
orderedPosLFC<-posLFC[order(posLFC$log2FoldChange, decreasing = TRUE), ]
padj_orderedPosLFC<-orderedPosLFC[order(orderedPosLFC$padj, decreasing = FALSE),]

# Find downregulated, and of those, significant genes
negLFC<-res[which(results_condition$log2FoldChange < 0),]
orderedNegLFC<-negLFC[order(negLFC$log2FoldChange, decreasing = FALSE), ]
padj_orderedNegLFC<-orderedNegLFC[order(orderedNegLFC$padj, decreasing = FALSE),]

# log normalize readcounts before multivariate analysis
rld <- rlogTransformation(dds, blind=FALSE)

# cluster heatmaps
library("pheatmap")
library("genefilter")
topGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)
mat <- assay(rld)[ topGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("sampleDepot", "sampleGenotype", "sizeFactor")])
pheatmap(mat, cluster_cols = TRUE, clustering_methods = "complete",annotation_col=df, fontsize = 10)
dev.copy(png, "deseq2_clusterHeatmap_test.png")
dev.off()

# Cooks Distance (to analyze outliers)
par(mar=c(15,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]),range=0, las=2)
dev.copy(png, "deseq2_cooksDistance_wt_RTvCE.png")
dev.off()

# Grab gene symbols from biomart
library(biomaRt)
gene_names <- row.names(res)
ensembl <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
genes.with.id=getBM(attributes=c("ensembl_gene_id", "external_gene_name"),values=gene_names, mart=ensembl)
annotationResults <- merge(genes.with.id, a s.data.frame(res), by.x="ensembl_gene_id", by.y="row.names")
write.csv(annotationResults, file = "Ing_WTvKO_20160115_annotated.csv")