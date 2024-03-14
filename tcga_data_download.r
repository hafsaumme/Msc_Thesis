library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")


query_TCGA = GDCquery(
  project = "TCGA-CHOL",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

GDCdownload(query = query_TCGA)

tcga_data_chol = GDCprepare(query_TCGA)
dim(tcga_data_chol)
assay(tcga_data_chol)
table(tcga_data_chol$sample_type)

saveRDS(tcga_data_chol,"E:/thesis/MSc/msc_figure/TCGA/tcga_data_chol.RDS")
