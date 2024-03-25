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
library("ggplot2")
library("ggpubr")


#======  Downloading TCGA LIHC data from TCGA database

query_TCGA = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

GDCdownload(query = query_TCGA)

tcga_data_lihc = GDCprepare(query_TCGA)

#====== Preparing TCGA dataframe with features===============

tcga_data_lihc$deceased <- ifelse(tcga_data_lihc$vital_status == "Alive", FALSE, TRUE)

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
tcga_data_lihc$overall_survival <- ifelse(tcga_data_lihc$vital_status == "Alive",
                                          tcga_data_lihc$days_to_last_follow_up,
                                          tcga_data_lihc$days_to_death)

tcga_data_lihc= tcga_data_lihc [,!is.na(tcga_data_lihc$overall_survival)]

#data normalization

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = assay(tcga_data_lihc),
                              colData = colData(tcga_data_lihc),
                              design = ~ 1)

# Removing genes with sum total of 10 reads across all samples

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# vst 
vsd <- vst(dds, blind=FALSE)

bulk.mtx=as.data.frame(assay(vsd)) %>%
  rownames_to_column(var = 'gene_id') %>%  
  left_join(., as.data.frame(tcga_data_lihc@rowRanges)[,c("gene_id","gene_name")], by = "gene_id") 


bulk.mtx$gene_id=NULL
rownames(bulk.mtx)=make.unique(bulk.mtx$gene_name)
bulk.mtx$gene_name=NULL
bulk.mtx_t=t(bulk.mtx)
bulk.mtx_t=rownames_to_column(as.data.frame(bulk.mtx_t),"sample_id")

clinical_data=as.data.frame(colData(tcga_data_lihc))
clinical_data=clinical_data[,c("sample_type","ajcc_pathologic_stage","gender","vital_status", "days_to_last_follow_up", "days_to_death","deceased","overall_survival")]

clinical_data$stage=clinical_data$ajcc_pathologic_stage
clinical_data$stage[which(str_detect(clinical_data$stage,"Stage III"))]="Stage III"
clinical_data$stage[which(str_detect(clinical_data$stage,"Stage IV"))]="Stage IV"
clinical_data$stage[which(is.na(clinical_data$stage))]="Normal"
clinical_data=rownames_to_column(clinical_data,"sample_id")
merge_lihc=merge(clinical_data, bulk.mtx_t,by="sample_id")



#survival curve MS4A1 gene

tumor_data=merge_lihc[(merge_lihc$sample_type == "Primary Tumor"),]

gene="MS4A1"

median_value <- median(tumor_data[, gene])

tumor_data[,gene] <- ifelse(tumor_data[,gene] >= median_value, "HIGH", "LOW")

tumor_data$Survival_years= tumor_data$overall_survival/365

fit <-survfit(as.formula(paste('Surv(Survival_years, deceased)~', gene)),data=tumor_data)

ggsurvplot(fit, data=tumor_data,surv.median.line = "none", pval=T, risk.table=F,
           xlab = "Time in years",legend.title="",legend = c(0.8, 0.9),legend.labs = 
             c("High expression", "Low expression"))$plot+ggtitle(gene)+
        theme(plot.title = element_text(hjust=0.5, face = "bold",size=13))


# Boxplot of MS4A1 gene between normal and tumor samples

ggboxplot(merge_lihc, 
          x = "sample_type", 
          y = "MS4A1", 
          fill="sample_type",
          outlier.shape = NA,
          bxp.errorbar = T,
          bxp.errorbar.width = 0.2,
          xlab=FALSE,ylab="log2(TPM + 1)",title = "MS4A1")+ 
  theme(legend.position = "none",plot.title = element_text(hjust=0.5, face = "bold",size=14))+
  geom_jitter(position = position_jitter(height = 4, width = .08),size=0.6)+
  geom_pwc(method = "t_test",label="p.format",bracket.nudge.y = 0.4, size=0.5,tip.length = 0.04)



# Boxplot for comparing stage of tumor

my_comparisons <- list( c("Stage I", "Stage II"),c("Stage I", "Stage III"), c("Stage I", "Stage IV") )


ggboxplot((tumor_data[which(!str_detect(tumor_data$stage,"Normal")),]), x = "stage", y = "MS4A1",
          fill = "stage")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)
