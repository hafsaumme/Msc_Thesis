# install devtools if necessary
BiocManager::install("TOAST")
devtools::install("TOAST")
# install the MuSiC package
devtools::install_github('xuranw/MuSiC')
# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("genefilter")
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
library(tidyverse)
library(SingleCellExperiment)
# load
library(MuSiC)
library(TOAST)
library(Biobase)
library(SingleCellExperiment)
library(MuSiC2)
BiocManager::install.packages(respos="https://CRAN.R-project.org/package=RefFreeEWAS")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
BiocManager::install("smplot2")
devtools::install_github('xuranw/MuSiC')
library(Seurat)
remove.packages("MuSiC")
#singlcell data



#bulk rna data

tcga_data <- readRDS("E:/thesis/MSc/msc_figure/TCGA/tcga_data.RDS")
clinical_data = colData(tcga_data)
library(DESeq2)
#start
dds <- DESeqDataSetFromMatrix(countData = assay(tcga_data),
                              colData = colData(tcga_data),
                              design = ~ 1)

# Removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


# vst 
vsd <- vst(dds, blind=FALSE)
brca_matrix_vst <- assay(vsd)
dim(assay(tcga_data))[1:10,1:10]
#end
group = factor(clinical_data$definition)
group = relevel(group, ref="Solid Tissue Normal")
design = model.matrix(~group)
dge = DGEList( # creating a DGEList object
  counts=assay(tcga_data),
  samples=colData(tcga_data),
  genes=as.data.frame(rowData(tcga_data)))

# filtering
keep = filterByExpr(dge,design) # defining which genes to keep
dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
rm(keep) 

dge = calcNormFactors(dge,method="TMM")

v = voom(dge,design,plot=TRUE)
View(dge$samples)
bulk.mtx = v$E
bulk.mtx = as.data.frame(assay(vsd)) %>%
  rownames_to_column(var = 'gene_id') %>%  
  left_join(., as.data.frame(tcga_data@rowRanges)[,c("gene_id","gene_name")], by = "gene_id") 
bulk.mtx["ALB",]
View(bulk.mtx)
intersect(bulk.mtx$gene_name,(obj.hcc_marker_100$X))
rownames(bulk.mtx)=make.names(bulk.mtx$gene_name,unique=TRUE)
bulk.mtx$gene_name=NULL
bulk.mtx$gene_id=NULL
#write.csv(bulk.mtx,"E:/thesis/MSc/msc_figure/TCGA/tcga_lihc_bulk.mtx.csv")
bulk.mtx=read.csv("E:/thesis/MSc/msc_figure/TCGA/tcga_lihc_bulk.mtx.csv",row.names = 1)
bulk.mtx_t=t(bulk.mtx)
View(as.data.frame(clinical_data))

write.csv(as.data.frame(clinical_data),"E:/thesis/MSc/msc_figure/TCGA/tcga_lihc_clinical_data_all.csv")
clinical_data=read.csv("E:/thesis/MSc/msc_figure/TCGA/tcga_lihc_clinical_data.csv",row.names = 1)


table(tcga_data$vital_status)
tcga_data$deceased <- ifelse(tcga_data$vital_status == "Alive", FALSE, TRUE)

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
tcga_data$overall_survival <- ifelse(tcga_data$vital_status == "Alive",
                                     tcga_data$days_to_last_follow_up,
                                     tcga_data$days_to_death)

View(as.data.frame(colData(tcga_data)))



clinical=as.data.frame(colData(tcga_data))[which(!is.na(tcga_data$overall_survival)),] 

which(colnames(clinical_brca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
clinical_data=clinical[,c("sample_type","gender","vital_status", "days_to_last_follow_up", "days_to_death","deceased","overall_survival")]
View(clinical_data)
View(CIBERSORTx_Job15_Results)

clinical_data$event <- ifelse(clinical_data$vital_status == "Alive", 0,1)
CIBERSORTx_Job15_Results <- (read.csv("E:/thesis/MSc/msc_figure/TCGA/cibersortx/CIBERSORTx_Job16_Results1.csv",row.names = 1))
CIBERSORTx_Job15_Results=CIBERSORTx_Job15_Results[,1:3]
gfg_standardized <- (CIBERSORTx_Job15_Results$B.cells.memory - mean(CIBERSORTx_Job15_Results$B.cells.memory)) / sd(CIBERSORTx_Job15_Results$B.cells.memory)
naive_mid <- median(CIBERSORTx_Job15_Results$Naive.Bcell)
naive <- ifelse(CIBERSORTx_Job15_Results$Naive.Bcell >= naive_mid, "HIGH", "LOW")

memory_mid <- median(CIBERSORTx_Job15_Results$Memory.Bcell)
memory <- ifelse(CIBERSORTx_Job15_Results$Memory.Bcell >= memory_mid, "HIGH", "LOW")
table(memory)
normal_index <- which(substr(colnames(rna),14,14) == '1')
tumor_index <- which(substr(colnames(rna),14,14) == '0')

# apply voom function from limma package to normalize the data
vm <- function(x){
  cond <- factor(ifelse(seq(1,dim(x)[2],1) %in% tumor_index, 1,  0))
  d <- model.matrix(~1+cond)
  x <- t(apply(x,1,as.numeric))
  ex <- voom(x,d,plot=F)
  return(ex$E)
}


clinical_data=rownames_to_column(clinical_data,"id")
CIBERSORTx_Job15_Results=rownames_to_column(CIBERSORTx_Job15_Results,"id")

CIBERSORTx_Job15_Results$Naive.Bcell=as.character(CIBERSORTx_Job15_Results$Naive.Bcell)
CIBERSORTx_Job15_Results$Memory.Bcell=as.character(CIBERSORTx_Job15_Results$Memory.Bcell)
CIBERSORTx_Job15_Results$Plasma.cell=as.character(CIBERSORTx_Job15_Results$Plasma.cell)
as.nu
clinical_data$Naive.Bcell=as.numeric(clinical_data$Naive.Bcell)
clinical_data$Memory.Bcell=as.numeric(clinical_data$Memory.Bcell)
clinical_data$Plasma.cell=as.numeric(clinical_data$Plasma.cell)

clinical_data$id <- gsub("-", ".", clinical_data$id)
clinical_data=merge(clinical_data,CIBERSORTx_Job15_Results,by="id")
View(clinical_data)
CIBERSORTx_Job15_Results$Naive.Bcell
typeof(clinical_data$overall_survival)
rownames(clinical_data)=clinical_data[,"id"]
clinical_data[,"id"]=NULL
res.cut <- surv_cutpoint(clinical_data, time = "overall_survival", event = "event",
                         variables = c("Memory.Bcell","Naive.Bcell","Plasma.cell"))

res.cat <- surv_categorize(res.cut)
clinical_data=rownames_to_column(clinical_data,"id")
res.cat =rownames_to_column(res.cat ,"id")
table(res.cat$Plasma.cell)

clinical_data[,"Naive_Bcell"]=clinical_data[,"Naive Bcell"]
clinical_data[,"Memory_Bcell"]=clinical_data[,"Memory Bcell"]
clinical_data[,"Plasma_cell"]=clinical_data[,"Plasma cell"]
write.csv(clinical_data,"E:/thesis/MSc/msc_figure/TCGA/cibersortx/clinicalData-_tcgaLICA")
res.cox <- coxph(Surv(overall_survival, deceased) ~ Memory.Bcell+Naive.Bcell+Plasma.cell,
                                                   data =  clinical_data)

Surv(overall_survival, deceased) ~ KPNA2 + CDCA8 + KIF2C + SFPQ + 
  KIF20A + NEIL3 + SOCS2 + CENPA + TRIP13 + UCK2

coeff=as.data.frame(res.cox$coefficients)
coeff=rownames_to_column(coeff,"subtype")
coeff[,"coefficients"]=coeff$`res.cox$coefficients`
coeff$`res.cox$coefficients`=NULL
  res.cox <- coxph(Surv(overall_survival, deceased) ~ (Memory_Bcell=="high")+(Naive_Bcell=="high") + 
                     (Plasma_cell=="high"), data =  clinical_data)
fit = survfit(Surv(overall_survival, deceased) ~ Memory.Bcell, data =  clinical_data)
ggsurvplot(fit, data=clinical_data)


library(ggpubr)

# Create a simple bar plot
ggbarplot(
  coeff, x = "subtype", y = "coefficients", 
  add = c("mean_sd", "jitter"),
  fill = "#BF504D"
)
devtools::install_github('smin95/smplot2')
library(tidyverse) 
library(cowplot) 
library(smplot2)
install.packages("smplot2")
ggplot(data =  coeff, mapping = aes(x = subtype, y = (coefficients), color = subtype)) +
  sm_bar() +
  scale_color_manual(values = sm_color('blue','orange',"green")) +
  ggtitle('A bar graph')
ggbarplot(
  coeff, x = "subtype", y = "coefficients",  
  add = c("mean_sd", "jitter"), 
  color = "supp", palette = c("#807F7F", "#BF504D","blue"),
  position = position_dodge(0.8)
)




boxplot(coefficients ~ subtype, data = coeff, col = "white")
stripchart(coefficients ~ subtype, data = coeff,
           method = "jitter",
           pch = 19,
           col = 2:4,
           vertical = TRUE,
           add = TRUE)
p=ggforest(res.cox,clinical_data )



data("ToothGrowth")
head(ToothGrowth)
ggerrorplot(ToothGrowth, x = "dose", y = "len", 
            desc_stat = "mean_sd")


CIBERSORTx_Job15_Results=read.csv("E:/thesis/MSc/msc_figure/TCGA/cibersortx/clinicalData-_tcgaLICA")
CIBERSORTx_Job15_Results=CIBERSORTx_Job15_Results[CIBERSORTx_Job15_Results$sample_type=="Primary Tumor" ,]
naive_mid <- median(CIBERSORTx_Job15_Results$Naive.Bcell.1)
CIBERSORTx_Job15_Results$Naive.Bcell.1 <- ifelse(CIBERSORTx_Job15_Results$Naive.Bcell.1 >= naive_mid, "HIGH", "LOW")

memory_mid <- median(CIBERSORTx_Job15_Results$Memory.Bcell)
memory <- ifelse(CIBERSORTx_Job15_Results$Memory.Bcell >= memory_mid, "HIGH", "LOW")



fit <-survfit(Surv(overall_survival, deceased)~Memory.Bcell.1,data=CIBERSORTx_Job15_Results)

fit



ggsurvplot(fit, data=CIBERSORTx_Job15_Results,surv.median.line = "none", pval=T, risk.table=F,
                    legend.title="",legend = c(0.8, 0.9),legend.labs = 
                      c("High expression", "Low expression"))$plot+ggtitle("Naive Bcells")
