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

meanSdPlot(assay(rld_sc))
meanSdPlot(assay(vsd_sc))
meanSdPlot(assay(ntd_sc))
meanSdPlot(assay(rld_dh))
meanSdPlot(assay(vsd_dh))
meanSdPlot(assay(ntd_dh))



