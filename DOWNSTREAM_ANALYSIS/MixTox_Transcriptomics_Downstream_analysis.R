#Mixture toxicity Project

#Downstream Transcriptomic analysis

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

#Install DESeq2 (version 1.30.1)
BiocManager::install("DESeq2")
library(DESeq2)

#IMPORT DATA

## STEP1 : Save the path of the HTSeq output folder in a variable

#for cerevisiae
Dir_sc <- "RAW_COUNT_DATA/HTSEQ_SC/"

#for hansenii
Dir_dh <- "RAW_COUNT_DATA/HTSEQ_DH/"

## STEP2 : Create a sample table with 3 columns. First column should be "sampleName" i.e. the name identifier of the samples. The second column should be "fileName" i.e. the file name list in this htseq output folder with the correct extentions, and the third column should be condition i.e. the treatment condition. The column names should be exactly same without the quotations. This is necessary for the functions to work properly. You can create this table in excel and then copy paste it in a .txt object. Later it can be read by the read.table function.

#for cerevisiae
sampleinfo_sc <- read.table("COMPILED_DATA/sampleTable_SC.txt", header = TRUE, sep = "\t", as.is = TRUE)

#for hansenii
sampleinfo_dh <- read.table("COMPILED_DATA/sampleTable_DH.txt", header = TRUE, sep = "\t", as.is = TRUE)

## STEP3 : Create the sampleTable object for both dataset

#for cerevisiae
sampleTable_sc <- data.frame(sampleName = sampleinfo_sc$sampleName,
                             fileName = sampleinfo_sc$fileName,
                             condition = sampleinfo_sc$condition)

#for hansenii
sampleTable_dh <- data.frame(sampleName = sampleinfo_dh$sampleName,
                             fileName = sampleinfo_dh$fileName,
                             condition = sampleinfo_dh$condition)

## STEP4 : Build the DESeqDataSet

#for cerevisiae
ddsHTSeq_sc <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_sc,
                                          directory = Dir_sc,
                                          design= ~ condition)

#for hansenii
ddsHTSeq_dh <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_dh,
                                          directory = Dir_dh,
                                          design= ~ condition)

## STEP5 : Print the object
#for cerevisiae
ddsHTSeq_sc
#for hansenii
ddsHTSeq_dh

##DATA TRANSFORMATION 

# STEP1: PRE-FILTERING
#for cerevisiae
keep<-rowSums(counts(ddsHTSeq_sc)) >= 10
dds_SC<-ddsHTSeq_sc[keep, ]
dds_SC_LOST<-ddsHTSeq_sc[!keep, ]

#for hansenii
keep<-rowSums(counts(ddsHTSeq_dh)) >= 10
dds_DH<-ddsHTSeq_dh[keep, ]
dds_DH_LOST<-ddsHTSeq_dh[!keep, ]

# STEP2: DATA TRANSFORMATION
#for cerevisiae
rld_sc <- rlogTransformation(dds_SC)
vsd_sc <- varianceStabilizingTransformation(dds_SC)
ntd_sc <- normTransform(dds_SC)
#for hansenii
rld_dh <- rlogTransformation(dds_DH)
vsd_dh <- varianceStabilizingTransformation(dds_DH)
ntd_dh <- normTransform(dds_DH)

# STEP3: DATA-VISUALIZATION > Effects of transformations on the variance

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("vsn")
library(vsn)

install.packages("cowplot")
library(cowplot)

pdf("R_FIGURES/1A_rld_sc.pdf")
meanSdPlot(assay(rld_sc))
dev.off()

pdf("R_FIGURES/1B_vsd_sc.pdf")
meanSdPlot(assay(vsd_sc))
dev.off()

pdf("R_FIGURES/1C_ntd_sc.pdf")
meanSdPlot(assay(ntd_sc))
dev.off()

pdf("R_FIGURES/2A_rld_dh.pdf")
meanSdPlot(assay(rld_dh))
dev.off()

pdf("R_FIGURES/2B_vsd_dh.pdf")
meanSdPlot(assay(vsd_dh))
dev.off()

pdf("R_FIGURES/2C_ntd_dh.pdf")
meanSdPlot(assay(ntd_dh))
dev.off()

#PCA PLOT
#for cerevisiae
pdf("R_FIGURES/3_PCA_SC.pdf", height = 4, width = 8)
plotPCA(rld_sc)
dev.off()

#for hansenii
pdf("R_FIGURES/4_PCA_DH.pdf", height = 4, width = 8)
plotPCA(rld_dh)
dev.off()

#DATA-VISUALIZATION 3 > SIMILARITY MATRIX

#Calculation of the euclidean distance of the transformed data. 
#We need to transpose the data since we want to calculate the distance between samples rather than among genes. 
#So we need the samples as rows rather than columns

library(DESeq2)
library(pheatmap)        
library(RColorBrewer)

#for cerevisiae
sampleDists_SC <- dist(t(assay(rld_sc)))
#Creating the corresponding matrix
sampleDistMatrix_SC <- as.matrix(sampleDists_SC)
#creating a distance matrix tagging the groups of the samples
sample_data_sc<-data.frame(Group=sampleinfo_sc$condition)
row.names(sample_data_sc)<-sampleinfo_sc$sampleName
#defining the colors we want to use depending on the condition
annotation_sc <- list(Group = c(Ctrl="green", Para="brown", Rapa="red", Salt="cyan", SaPa="blue", SaRa="magenta"))

#for hansenii
sampleDists_DH <- dist(t(assay(rld_dh)))
#Creating the corresponding matrix
sampleDistMatrix_DH <- as.matrix(sampleDists_DH)
#creating a distance matrix tagging the groups of the samples
sample_data_dh<-data.frame(Group=sampleinfo_dh$condition)
row.names(sample_data_dh)<-sampleinfo_dh$sampleName
#defining the colors we want to use depending on the condition
annotation_dh <- list(Group = c(Ctrl="green", Para="brown", Rapa="red", Salt1="cyan", Salt2="turquoise4", SaPa="blue", SaRa="magenta"))

#Make a color palette for heatmap
colfunc<-colorRampPalette(c("goldenrod4", "goldenrod", "white", "turquoise", "turquoise4"))

#Plotting the heatmap for cerevisiae 
pdf("R_FIGURES/5_Similarity_mat_SC.pdf")
pheatmap(sampleDistMatrix_SC, 
         annotation_col = sample_data_sc, 
         annotation_colors = annotation_sc, 
         color = colfunc(50))
dev.off()


par(mfcol=c(4,2))
hist(salt_SC$log2FoldChange,
     breaks = 50,
     xlab = "log2FoldChange", 
     ylab = "Frequency",
     main = "Control vs Salt(600mM)",
     col = "skyblue",
     xlim = c(-6, 6),
     ylim = c(0, 3000))
plot(NULL, 
     xlim=c(-6,6), 
     ylim=c(0,3000), 
     ylab="", 
     xlab="")
hist(para_SC$log2FoldChange,
     breaks = 50,
     xlab = "log2FoldChange", 
     ylab = "Frequency",
     main = "Control vs Paraquat (1000 ug/ml)",
     col = "skyblue",
     xlim = c(-6, 6),
     ylim = c(0, 3000))
hist(rapa_SC$log2FoldChange,
     breaks = 50,
     xlab = "log2FoldChange", 
     ylab = "Frequency",
     main = "Control vs Rapamycin (1 ug/ml)",
     col = "skyblue",
     xlim = c(-6, 6),
     ylim = c(0, 3000))
par(mfcol=c(4,2))
hist(salt1_DH$log2FoldChange,
     breaks = 50,
     xlab = "log2FoldChange", 
     ylab = "Frequency",
     main = "Control vs Salt(650mM)",
     col = "skyblue",
     xlim = c(-6, 6),
     ylim = c(0, 3000))
hist(salt2_DH$log2FoldChange,
     breaks = 50,
     xlab = "log2FoldChange", 
     ylab = "Frequency",
     main = "Control vs Salt(1200mM)",
     col = "skyblue",
     xlim = c(-6, 6),
     ylim = c(0, 3000))
hist(para_DH$log2FoldChange,
     breaks = 50,
     xlab = "log2FoldChange", 
     ylab = "Frequency",
     main = "Control vs Paraquat (500 ug/ml)",
     col = "skyblue",
     xlim = c(-6, 6),
     ylim = c(0, 3000))
hist(rapa_DH$log2FoldChange,
     breaks = 50,
     xlab = "log2FoldChange", 
     ylab = "Frequency",
     main = "Control vs Rapamycin (1 ug/ml)",
     col = "skyblue",
     xlim = c(-6, 6),
     ylim = c(0, 3000))
