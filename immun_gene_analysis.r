library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(tibble)
library("stringr")
library(stringr)
library(glmGamPoi)
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
obj.hcc_icca <- readRDS("E:/thesis/MSc/msc_figure/hcc_icca/obj.hcc_icca.rds")
DefaultAssay(obj.hcc_icca)="RNA"
obj.hcc=subset(obj.hcc_icca, subset=orig.ident=="HCC")
VlnPlot(object = obj.hcc, 
        features = c("C3","GC","CPB2","FGB","HRG","KNG1","SERPINC1"),ncol2)
#ggsave("E:/thesis/MSc/msc_figure/TCGA/hcc_hubgene_vlnplot.jpg", width = 12, height = 13, units = c("in"), dpi = 300)

type=obj.hcc$Type

rownames(res_chol)
obj.hcc <- CreateSeuratObject(counts = obj.hcc@assays$RNA@counts, 
                              project = "hcc")


obj.hcc=SetIdent(obj.hcc, value=type)

obj.hcc <- NormalizeData(obj.hcc)
obj.hcc  <- FindVariableFeatures(obj.hcc, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(obj.hcc)
obj.hcc <- ScaleData(obj.hcc, features = all.genes)


obj.hcc_marker <- FindMarkers(obj.hcc, ident.1 = "Malignant cell",only.pos=T,test.use = "wilcox", logfc.threshold = 0.25)
#write.csv(obj.hcc_marker,"E:/thesis/MSc/msc_figure/hcc_icca/obj.hcc_marker.csv")
obj.hcc_marker=read.csv("E:/thesis/MSc/msc_figure/hcc_icca/obj.hcc_marker.csv")

img=read.csv("E:/thesis/MSc/msc_figure/TCGA/InnateDB_genes.csv")
img_gene=img$name



#tcga
tcga_data <- readRDS("E:/thesis/MSc/msc_figure/TCGA/tcga_data.RDS")
tcga_data_tumot=tcga_data [,colData(tcga_data)$sample_type == "Primary Tumor"]
dim(tcga_data)
tcga_data_tumot=tcga_data
tcga_data_tumot$deceased <- ifelse(tcga_data_tumot$vital_status == "Alive", FALSE, TRUE)

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
tcga_data_tumot$overall_survival <- ifelse(tcga_data_tumot$vital_status == "Alive",
                                           tcga_data_tumot$days_to_last_follow_up,
                                           tcga_data_tumot$days_to_death)

tcga_data_tumot= tcga_data_tumot [,!is.na(tcga_data_tumot$overall_survival)]

group = factor(colData(tcga_data_tumot)$definition)
group = relevel(group, ref="Solid Tissue Normal")
design = model.matrix(~group)
dge = DGEList( # creating a DGEList object
  counts=assay(tcga_data_tumot),
  samples=colData(tcga_data_tumot),
  genes=as.data.frame(rowData(tcga_data_tumot)))

# filtering
keep = filterByExpr(dge,design) # defining which genes to keep
dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
rm(keep) 
dge$genes$length=nchar(dge$genes$gene_id)
#dge = calcNormFactors(dge,method="TMM")
#========tpm norm
tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

g <- data.frame( ensembl_gene_id = dge$genes$gene_id , 
                 transcript_length = dge$genes$length,
                 stringsAsFactors = F, dge$genes$gene_id)
g <- g[!duplicated(g$ensembl_gene_id),]
rownames(g)=g$ensembl_gene_id
cf=dge$counts
igenes <- intersect(rownames(cf),g$ensembl_gene_id)
g1 <- g[igenes,]
cf1 <- cf[igenes,]
all.equal(rownames(cf1),g1$ensembl_gene_id)
library(DESeq2)
ct <- tpm(cf1,g1$transcript_length)
ct <- log2( ct + 1 )
boxplot(ct,ylab=expression('Log'[2]~'Read counts'),las=2,main="TPM")
#end
tpm=tmm(dge)
v=log(counts(dge$counts, normalized = TRUE)+1)
logcounts <- log2( counts(dds, normalized=TRUE) + 1 )
View(dge)
v = voom(dge,design,plot=TRUE)
View(as.data.frame(dge$samples$norm.factors))

data=v$E
tcga_data_tumot= tcga_data_tumot [,colData(tcga_data_tumot)$sample_type == "Primary Tumor"]

bulk.mtx=as.data.frame(data) %>%
  rownames_to_column(var = 'gene_id') %>%  
  left_join(., as.data.frame(tcga_data_tumot@rowRanges)[,c("gene_id","gene_name")], by = "gene_id") 


bulk.mtx$gene_id=NULL

rownames(bulk.mtx)=make.unique(bulk.mtx$gene_name)
bulk.mtx$gene_name=NULL
bulk.mtx_t=t(bulk.mtx)
bulk.mtx_t=rownames_to_column(as.data.frame(bulk.mtx_t),"sample_id")
clinical_data=as.data.frame(colData(tcga_data_tumot))
clinical_data=clinical_data[,c("sample_type","gender","vital_status","tissue_or_organ_of_origin", "days_to_last_follow_up", "days_to_death","deceased","overall_survival")]
clinical_data=rownames_to_column(clinical_data,"sample_id")

merge_lihc=merge(clinical_data, bulk.mtx_t,by="sample_id")
dim(merge_lihc)
#View(bulk.mtx)
#write.csv(common,"E:/thesis/MSc/msc_figure/hcc_icca/immune_deg_common.csv")
#merge_lihc=read.csv("E:/thesis/MSc/msc_figure/hcc_icca/merge_icc_dds.csv")
merge_lihc=read.csv("E:/thesis/MSc/msc_figure/TCGA/tcga_hcc_tpm.csv")

obj.hcc_marker=read.csv("E:/thesis/MSc/msc_figure/hcc_icca/obj.hcc_marker.csv")
merge_lihc=read.csv("E:/thesis/MSc/msc_figure/TCGA/tcga_hcc_tpm.csv")
merge_lihc=merge_lihc[(merge_lihc$sample_type == "Primary Tumor"),]

common=intersect(img_gene, obj.hcc_marker$X)
common=intersect(common,colnames(merge_lihc))
common=common[str_detect(common,'MT-') == FALSE]

covariates= rownames(res_chol)
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(overall_survival, deceased)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = merge_lihc)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res_chol <- t(as.data.frame(univ_results, check.names = FALSE))

res_chol=filter(as.data.frame(res_chol), p.value<0.05)
Cluste1MCC_top38_default_node <- read_csv("E:/thesis/MSc/msc_figure/cytoscape/Cluste1MCC_top38 default node.csv")
hubgene <- read_csv("E:/thesis/MSc/msc_figure/cytoscape/hubgene22.csv")
hubgene=hubgene$`query term`
cytoscape_gene=Cluste1MCC_top38_default_node$`query term`
res_chol=filter(as.data.frame(res_chol), rownames(res_chol)%in%hubgene)
#write.csv(res_chol,"E:/thesis/MSc/msc_figure/cytoscape/cytoscapeGeneImmunehub.csv")
rownames(res_chol)
res_chol=read.csv("E:/thesis/MSc/msc_figure/cytoscape/cytoscapeGeneImmunehub.csv",row.names = 1)
#=====AIC============
options(scipen = 999)
library(glmnet)
common_pos=paste(rownames(res_chol), collapse="+")
common_pos
X <- model.matrix(overall_survival ~  C3+GC+CPB2+FGB+FGA+HRG+KNG1+SERPINC1, data=merge_lihc)[,-1]
#and the outcome
Y <- merge_lihc[,"overall_survival"] 

cv.lambda.lasso <- cv.glmnet(x=X, y=Y, 
                             alpha = 1,nfolds = 10,type.measure = "deviance") 
plot(cv.lambda.lasso)

l.lasso.min <- cv.lambda.lasso$lambda.min
lasso.model <- glmnet(x=X, y=Y,
                      alpha  = 1, 
                      lambda = l.lasso.min)
lasso.model$beta 
lasso.model$
summary(lasso.model)
ols.model <- glm(overall_survival~  CPB2+HRG+SERPINC1, data=merge_lihc)

summary(ols.model)     
#end=======

library(ggplot2)
#install.packages("forestmodel")
common_pos=paste((c("C3","GC","CPB2","FGB","HRG","KNG1","SERPINC1")), collapse="+")
common_pos= "HRG+CPB2+SERPINC1"
f1=as.formula(paste("Surv(overall_survival, deceased) ~ ",
                    common_pos))

gene.cox <- coxph(f1,data =  (merge_lihc) )
sm_cox=summary(gene.cox)
sm_cox
gene.cox$
#suvf=survfit(f1,data =  (merge_lihc) )

library(forestmodel)
panels <- list(
  list(width = 0.03),
  list(width = 0.09,display = ~variable, fontface = "bold", heading = "Gene"),
  
  
  
  list(width = 0.02, item = "vline", hjust = 0.5),
  list(
    width = 1, item = "forest", hjust = 0.5, heading = "Hazard ratio", linetype = "dashed",
    line_x = 0
  ),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.3,heading ="Reference", display = ~ ifelse(reference, "Reference", sprintf(
    "%0.2f (%0.2f, %0.2f)",
    trans(estimate), trans(conf.low), trans(conf.high)
  )), display_na = NA)
  ,list(width = 0.015, item = "vline", hjust = 0.5),
  list(
    width = 0.02,
    display = ~ ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
    display_na = NA, hjust = 1, heading = "p"
  ),
  list(width = 0.015)
)

library(tidymodels)
library(dplyr)
library("survival")
#install.packages("tidymodels")
res_chol=as.double(res_chol)
forestmodel::forest_model(model_list=univ_models,covariates=rownames(res_chol),panels,merge_models =TRUE)


#ggsave("E:/thesis/MSc/msc_figure/TCGA/hcc_forest_univariate.jpg", width = 10, height = 3.5, units = c("in"), dpi = 300)

#===========================================================
# get median value
library(dplyr)
tcga_data <- readRDS("E:/thesis/MSc/msc_figure/TCGA/tcga_data.RDS")
clinical_data = colData(tcga_data)
View(as.data.frame(clinical_data))
clinical_data$stage=clinical_data$ajcc_pathologic_stage
clinical_data=select(as.data.frame(clinical_data), stage)
clinical_data$stage[which(str_detect(clinical_data$stage,"Stage III"))]="Stage III"
clinical_data$stage[which(str_detect(clinical_data$stage,"Stage IV"))]="Stage IV"
clinical_data=rownames_to_column(clinical_data,"sample_id")
clinical_data$stage[which(is.na(clinical_data$stage))]="Normal"
#clinical_data= clinical_data[which(!is.na(clinical_data$stage)),]


clinical_data=read.csv("E:/thesis/MSc/msc_figure/TCGA/clinical_data_lihc_stage.csv")
MS4A1=select(readRDS("E:/thesis/MSc/msc_figure/TCGA/tcga_hcc_vsd.rds"),sample_id,MS4A1,CD79A,TCL1A)
MS4A1$Naive_Bcell=MS4A1$MS4A1+MS4A1$CD79A+MS4A1$TCL1A
MS4A1=left_join(clinical_data,MS4A1, by="sample_id")
immunecell=(read.csv("E:/thesis/MSc/msc_figure/TCGA/cibersortx/CIBERSORTx_Job34_Results.csv",row.names = 1))
immunecell=rownames_to_column(immunecell,"sample_id")
immunecell$sample_id <- gsub("\\.", "-", immunecell$sample_id)

merge_lihc=readRDS("E:/thesis/MSc/msc_figure/TCGA/tcga_hcc_vsd.rds")
merge_lihc=filter(merge_lihc, sample_id%in% sig_id)
merge_lihc=read.csv("E:/thesis/MSc/msc_figure/TCGA/tcga_hcc_tpm.csv")
merge_lihc=read.csv("E:/thesis/MSc/msc_figure/TCGA/tcga_hcc_dge.csv")
merge_lihc <-(read.csv("E:/thesis/MSc/msc_figure/TCGA/cibersortx/clinicalData-_tcgaLICA",row.names = 1))
merge_lihc$id <- gsub("\\.", "-", merge_lihc$id)
merge_lihc$sample_id=merge_lihc$id

merge_lihc1=merge(immunecell,merge_lihc, by="sample_id")
cibersort <-(read.csv("E:/thesis/MSc/msc_figure/TCGA/cibersortnaiv.csv"))
cibersort$sample_id <- gsub("\\.", "-", cibersort$sample_id)
merge_lihc1=merge(merge_lihc1,merge_lihc2, by="sample_id")
merge_lihc1=filter(merge_lihc1, P.value<0.05)
merge_lihc1$B.cells.naive

merge_lihc$sample_type[which(str_detect(merge_lihc$sample_type,"Tumor"))]="Tumor"
merge_lihc$Naive_Bcell=merge_lihc$Naive_Bcell/3
merge_lihc$ms=arrange(merge_lihc, stage)
merge_lihc$Sample=merge_lihc$sample_type
merge_lihc=merge_lihc[(merge_lihc$sample_type == "Primary Tumor"),]
sig_id=merge_lihc$sample_id
merge_lihc2=merge_lihc
ggbarplot(
  merge_lihc, x = "sample_type", y = "FGB",
  add = c("mean_sd"),
  fill = c("#807F7F"),ylab="Cox coefficient",xlab="",width = 0.4
)+ stat_compare_means(method = "t.test")geom_jitter(position = position_jitter(height = .1, width = .08),size=1)+
  c("Stage I", "Stage II")
table(merge_lihc$sample_type)
 
 my_comparisons <- list(    c("Stage I", "Stage II"),c("Stage I", "Stage III"), c("Stage I", "Stage IV") )


 
 ggboxplot(merge_lihc, 
          x = "sample_type", 
          y = "CPB2", 
         
           fill="sample_type",
          outlier.shape = NA,
          bxp.errorbar = T,
          bxp.errorbar.width = 0.2,
          
          
          
          
          xlab=FALSE,ylab="log2(TPM + 1)",title = "CPB2")+ theme(legend.position = "none",plot.title = element_text(hjust=0.5, face = "bold",size=14))+
  geom_jitter(position = position_jitter(height = 4, width = .08),size=0.6)+geom_pwc(method = "t_test",label="p.format",bracket.nudge.y = 0.4, size=0.5,tip.length = 0.04)

  stat_compare_means(label = "p.format",   bracket.size=0.5, label.y = c(17),label.x = 1.3)

 tumor=merge_lihc$C3[merge_lihc$sample_type=="Tumor"]
 normal=merge_lihc$C3[merge_lihc$sample_type=="Normal"]
 t.test(tumor, normal)
 ggsave("E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/trajectory/Final/CPB2_sample.jpg", width = 4, height = 4.5, units = c("in"), dpi = 300)
 


merge_lihc[, "MS4A1"]
merge_lihc=merge_lihc1
gene="NaiveB"
#merge_lihc= merge_lihc[which(!is.na(merge_lihc$sample_type)),]
median_value <- median(merge_lihc[, gene])

# denote which cases have higher or lower expression than medain count
merge_lihc[,gene] <- ifelse(merge_lihc[,gene] >= median_value, "HIGH", "LOW")

#coxph(Surv(overall_survival, deceased)~C1S, merge_lihc) 

# fitting survival curve -----------
merge_lihc$Survival_years= merge_lihc$overall_survival/365
#survfit(as.formula(paste('Surv(Survival_years, deceased)~', gene)),data=merge_lihc)


fit <-survfit(as.formula(paste('Surv(Survival_years, deceased)~', gene)),data=merge_lihc)

fit



ggsurvplot(fit, data=merge_lihc,surv.median.line = "none", pval=T, risk.table=F,
                  xlab = "Time in years",legend.title="",legend = c(0.8, 0.9),legend.labs = 
                    c("High expression", "Low expression"))$plot+ggtitle(gene)+
  theme(plot.title = element_text(hjust=0.5, face = "bold",size=13))



C3 GC CPB2 FGB, FGA FGG KLKB1 F12 HRG KNG1 PROC PLG SERPINC1
, GC, CPB2, FGB, FGA, FGG, KLKB1, F12, HRG, KNG1, PROC, PLG, SERPINC1

library(ggpubr)
library(ggplot2)
ggarrange(C3, GC, CPB2, FGB, FGA, FGG, KLKB1, F12, HRG, KNG1, PROC, PLG, SERPINC1,nrow=5,ncol=3)


ggsave("E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/trajectory/Final/riskgroup.jpg", width = 5.5, height = 4.5, units = c("in"), dpi = 300)


#correlation analysis


immunecell=(read.csv("E:/thesis/MSc/msc_figure/TCGA/cibersortx/CIBERSORTx_Job16_Results1.csv",row.names = 1))
merge_lihc=read.csv("E:/thesis/MSc/msc_figure/TCGA/tcga_hcc_tpm.csv")
merge_lihc_hub=merge_lihc[,c("sample_id","C3","GC","CPB2","FGB","HRG","KNG1","SERPINC1")]
library(stringr)
immunecell=rownames_to_column(immunecell,"sample_id")
immunecell$sample_id=gsub("\\.", "-", immunecell$sample_id)
merge_immunecell=merge(immunecell,merge_lihc_hub,by="sample_id")
merge_immunecell=merge_immunecell%>%select(-c("P.value","Correlation","RMSE","sample_id"))
merge_lihc$sample_id

result = cor(merge_immunecell,method = "spearman")
result =result[,1:22]
result=result[23:29,]
result=round(result,2)


install.packages("reshape2")
library(reshape2)
melted_corr_mat <- melt(result)

ggplot(data = melted_corr_mat, aes(x=Var1, y=Var2, 
                                   fill=value)) + 
  geom_tile() +
  geom_text(aes(Var1, Var2, label = value), 
            color = "black", size = 3)+
scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0,  space = "Lab",
                       name="Spearman\ncorrelation")+
  theme(axis.text.x = element_text(angle=60,hjust=1,size=11,color="black"),
        axis.text.y = element_text(size=11,color="black"),
        axis.title = element_blank()
        )+coord_flip()


#ggsave("E:/thesis/MSc/msc_figure/TCGA/hcc_hubgeneSelected_correlation.jpg", width = 10, height = 5, units = c("in"), dpi = 300)

sm_cox_df=as.data.frame(sm_cox$coefficients)
sm_cox_df_coeff=as.data.frame(t(select(sm_cox_df,coef)))
rownames(sm_cox_df_coeff)=NULL


#===============Calculating risk score from coxph significant gene==========

merge_lihc$riskScore=(merge_lihc$HRG*sm_cox_df_coeff$HRG)+(merge_lihc$CPB2*sm_cox_df_coeff$CPB2)+(merge_lihc$SERPINC1*sm_cox_df_coeff$SERPINC1)

merge_lihc$Survival_years= merge_lihc$overall_survival/365
median_value_rs = median(merge_lihc$riskScore)
merge_lihc$RS_Group = ifelse(  merge_lihc$riskScore >= median_value_rs,  "High risk","Low risk")
View(select(merge_lihc,riskScore,RS_Group))
summary(fit_rs)
table(merge_lihc$RS_Group)

fit_rs = survfit(Surv(Survival_years, deceased) ~ RS_Group, data=merge_lihc)
ggsurv=ggsurvplot(fit_rs, data=merge_lihc,surv.median.line = "none", pval=F, risk.table=F,
                  xlab = "Time in years",legend.labs = 
                    c("High risk group", "Low risk group"),legend.title="")


ggsurv$plot+ ggplot2::annotate("text", 
                               x = 3, y = 0.09, # x and y coordinates of the text
                               label = "High risk group (median) - 4.44 years\nLow risk group (median) - 5.80 years", size = 3.5)
merge_lihc$lp=predict(gene.cox,type="lp")
library(survivalROC)
library(timeROC)
model=Surv(merge_lihc$Survival_years, merge_lihc$deceased)
ROC.bili.cox<-timeROC(T=merge_lihc$Survival_years,
                      delta=merge_lihc$deceased,
                      marker=-merge_lihc$RS_Group,
                      cause=1,weighting="cox",
                      times=c(1,2,3,4,5,6,7,8,9,10))
plot(ROC.bili.cox, time=6)

Mayo4.2= survivalROC(Stime=merge_lihc$Survival_years,  
                     status=merge_lihc$deceased,      
                     marker = merge_lihc$riskScore,     
                     predict.time = 10 , method="KM")
plot(Mayo4.2$FP, Mayo4.2$TP, type="l") 
