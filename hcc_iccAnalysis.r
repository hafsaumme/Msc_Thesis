
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(tibble)
library("stringr")
library(stringr)
library(glmGamPoi)


obj.hcc_icca <- readRDS("E:/thesis/MSc/msc_figure/hcc_icca/obj.hcc_icca.rds")
DefaultAssay(obj.hcc_icca)="RNA"
obj.hcc=subset(obj.hcc_icca, subset=orig.ident=="HCC")
type=obj.hcc$Type

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
obj.hcc_marker_sig  = obj.hcc_marker  %>%arrange(p_val_adj)%>%top_n(200)
 
obj.hcc_marker_200  = obj.hcc_marker  %>%
  slice_max(n = 300, order_by = (avg_log2FC))
obj.hcc_marker=filter(obj.hcc_marker,p_val_adj<=0.01,avg_log2FC>=0.5)
obj.hcc_icca=SetIdent(obj.hcc_icca, value= obj.hcc_icca$Type)
obj.hcc_icca <- FindMarkers(obj.hcc_icca, ident.1 = "Malignant cell",only.pos=T,test.use = "wilcox", logfc.threshold = 0.25)
obj.hcc_icca  = obj.hcc_icca  %>%
  slice_max(n = 50, order_by = (avg_log2FC))

malignant=subset(obj.hcc_icca, subset=Type=="Malignant cell")
table(malignant$orig.ident)

malignant_marker=FindAllMarkers(malignant, only.pos=T,test.use = "wilcox", logfc.threshold = 0.25)
malignant_marker=read.csv("E:/thesis/MSc/msc_figure/hcc_icca/malignant_marker_hcc_icc.csv")
malignant_marker = filter(malignant_marker, p_val_adj<=0.01)
malignant_marker = malignant_marker %>%
  group_by(cluster) %>%
  slice_max(n = 200, order_by = (avg_log2FC))



hcc_malignant_marker=malignant_marker[malignant_marker$cluster=="GSE125449.set1_2_HCC",]
icca_malignant_marker=malignant_marker[malignant_marker$cluster=="GSE125449.set1_2_iCCA",]

hcc_malignant_marker <- FindMarkers(malignant, ident.1 = "GSE125449.set1_2_HCC",only.pos=T,test.use = "wilcox", logfc.threshold = 0.25)

library(dplyr)
library(stringr)
library(ggplot2)
library(fgsea)
library(GO.db)
library(ggfittext)
library(Seurat)
library(clusterProfiler)
library(purrr)
library(data.table)
require(DOSE)
library(DOSE)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
BiocManager::install("clusterProfiler.rdb")
#

#BiocManager::install("gage")
library(org.Hs.eg.db)
SYMBOL2EG <-
  eval(parse(text = sprintf(
    'org.%s.egSYMBOL2EG', 'Hs'
  )))
logFC_score <- obj.hcc_marker$avg_log2FC
logFC_score
de_genes<- (obj.hcc_marker$X)
length(x = de_genes)
names(logFC_score) <- de_genes
logFC_score

genes <- intersect(de_genes, mappedkeys(SYMBOL2EG)) #  Get the gene symbol that are mapped to an entrez gene identifiers
#length(x = genes)


logFC_score<- logFC_score[genes] # access logFC score 
#length(x = logFC_score) 
gene_entrez <-genes %>% SYMBOL2EG[.] %>% as.list %>% map( ~ .[1]) %>% purrr::simplify()
#length(x = gene_entrez)

names(logFC_score) <- gene_entrez
logFC_score # input for rank in fgsea


library(data.table)
require(DOSE)
library(DOSE)
library(clusterProfiler)
library(enrichplot)
library(ggupset)

BiocManager::install("clusterProfiler")
biocLite("clusterProfiler", siteRepos = "https://bioconductor.org/packages/release/bioc/")
gse <- gseGO(geneList=sort(logFC_score, decreasing = T), 
             ont ="BP", 
             scoreType = "pos",
            
            
             OrgDb = org.Hs.eg.db
)

View(gse1@result)
#remove.packages("clusterProfiler")

#saveRDS(gse,"E:/thesis/MSc/msc_figure/hcc_icca/gse_hcc_malignant.rds")
gse=readRDS("E:/thesis/MSc/msc_figure/hcc_icca/gse_hcc_malignant.rds")
gse1=filter(gse, !str_detect(Description,"positive|tube|blood|protein|plasma"))
gse1=(gse1%>% arrange(desc(enrichmentScore)))
gse1_20=gse1[1:20]

dotplot(gse1, showCategory=20) 
barplot(gse1, showCategory=20) 
barplot((gse1), showCategory=20) 
dotplot(gse, x="enrichmentScore", showCategory=20)

  
  

  gse1_20=read.csv("E:/thesis/MSc/msc_figure/hcc_icca/gsego_hcc20.csv")
  
  gse1_20$GeneRatio=(round(-log10(gse1_20$enrichmentScore),2))
  
  gse1_20$GeneRatio= gse1_20$enrichmentScore
  gse1_20=gse1_20%>%arrange(GeneRatio,decreasing=T)
 typeof(gse1_20$GeneRatio)
  ggplot(gse1_20, 
         aes(x=reorder(Description,enrichmentScore), 
             y= enrichmentScore,
             fill="Upregulated")) +  labs(title = paste("GO"))+
    
    geom_bar(stat = "identity",width = 0.6, show.legend = FALSE)+coord_flip()+
    theme(panel.background = element_blank(),
          axis.line.x.top = element_line(size = 1),
          axis.line.x.bottom = element_line(size = 0.3),
          axis.line.y.left = element_line(linetype="dotted"),
          axis.text=element_text(color="black",size=10),
          axis.title.y =element_blank(),
          plot.title = element_text(color="black",size=10,face = "bold",hjust = 0.5) )

  gsekegg1_20$p.adjust
  
  ggplot(gsekegg1_20, 
         aes(x=reorder(Description,GeneRatio), 
             y= GeneRatio,
    fill = p.adjust
  )) +
    geom_col(width = 0.8)+coord_flip()+ 
    scale_fill_gradient(low = "red", high = "blue")+theme_bw()+
    
    scale_x_discrete(labels = function(x) str_wrap(x, width =70))+
    theme(axis.text=element_text(color="black",size=12),
          panel.border = element_rect(color="black",fill=NA,linewidth =0.5),
    )+
    labs(title = "",
         x="",
         y = 'GeneRatio')
  #ggsave("E:/thesis/MSc/msc_figure/hcc_icca/hcc_malignantcell_kegg_enrichment_final.jpg", width = 7, height = 4.5, units = c("in"), dpi = 300)
  
  
  
 gsekegg <- gseKEGG(geneList     = sort(logFC_score, decreasing = T),
                 organism     = 'hsa',
                 scoreType = "pos",
                 
                 minGSSize = 1,
                 maxGSSize = Inf,
                 pvalueCutoff = 5,
                
                 verbose      = FALSE)
 View(gsekegg1@result)
 
#saveRDS(gsekegg,"E:/thesis/MSc/msc_figure/hcc_icca/gsekegg_hcc_malignant.rds")
 gsekegg=readRDS("E:/thesis/MSc/msc_figure/hcc_icca/gsekegg_hcc_malignant.rds")
 gsekegg1=filter(gsekegg, !str_detect(Description,"Coronavirus|protein"),pvalue<0.02)
 gsekegg1=filter(gsekegg1,ID %in% c("hsa04610","hsa04979","hsa04918","hsa00980","hsa05150",
 "hsa05143","hsa04978","hsa05204","hsa05204","hsa05208","hsa03320","hsa04260","hsa04611"))
 gsekegg1=(gsekegg1%>% arrange(desc(enrichmentScore)))
 gsekegg1_20=gsekegg1[1:12]
 #options(scipen = 999)
 gsekegg1_20$p.adjust=gsekegg1_20$pvalue
 gsekegg1_20=gsekegg1_20%>%arrange(GeneRatio,decreasing=T)
 gsekegg1_20$GeneRatio= gsekegg1_20$enrichmentScore
 typeof(gsekegg1_20$GeneRatio)
 ggplot(gsekegg1_20, 
        aes(x=reorder(Description,GeneRatio), 
            y=GeneRatio,
            fill="Upregulated")) +  labs(title = paste("KEGG"))+
   geom_bar(stat = "identity",width = 0.6, show.legend = FALSE)+coord_flip()+
   theme(panel.background = element_blank(),
         axis.line.x.top = element_line(size = 1),
         axis.line.x.bottom = element_line(size = 0.3),
         axis.line.y.left = element_line(linetype="dotted"),
         axis.text=element_text(color="black",size=10),
         axis.title.y =element_blank(),
         plot.title = element_text(color="black",size=10,face = "bold",hjust = 0.5) )
 
 
 
 #==============Survival analysis==========================================================================
 #ggsave("E:/thesis/MSc/msc_figure/hcc_icca/hcc_malignantcell_KEGG_enrichment.jpg", width = 6, height = 4, units = c("in"), dpi = 300)
 LIHC_mrna <- read_excel("E:/thesis/MSc/msc_figure/TCGA/cibersortx/LIHC_mrna.xlsx")
 genes=LIHC_mrna$`TCGA Name`
 
 ##cox coefficent of hcc gene
 #bulk.mtx=read.csv("E:/thesis/MSc/msc_figure/TCGA/tcga_lihc_bulk.mtx.csv",row.names = 1)
 #bulk.mtx_t=t(bulk.mtx)
 #bulk.mtx_t=rownames_to_column(as.data.frame(bulk.mtx_t),"sample_id")
 
 #clinical_data=read.csv("E:/thesis/MSc/msc_figure/TCGA/tcga_lihc_clinical_data.csv",row.names = 1)
 #clinical_data=rownames_to_column(clinical_data,"sample_id")
 #clinical_data$sample_id= gsub('\\-', '.', clinical_data$sample_id)
 #clinical_data_type=select(clinical_data,deceased, sample_id,sample_type)
 
 #bulk.mtx_t=merge(clinical_data_type, bulk.mtx_t,by="sample_id")
 #bulk.mtx_t=bulk.mtx_t[bulk.mtx_t$sample_type=="Primary Tumor",]
 #bulk.mtx_t=bulk.mtx_t[bulk.mtx_t$deceased==TRUE,]
 #rownames(bulk.mtx_t)=bulk.mtx_t$sample_id
 #bulk.mtx_t$sample_id=NULL
 #bulk.mtx_t$sample_type=NULL
 #bulk.mtx_t$deceased=NULL
 #dim(bulk.mtx_t)
 #rownames(bulk.mtx_t)= gsub('\\.', '-', rownames(bulk.mtx_t))
 
 
 
 tcga_data <- readRDS("E:/thesis/MSc/msc_figure/TCGA/tcga_data.RDS")
 tcga_data_tumot=tcga_data [,colData(tcga_data)$sample_type == "Primary Tumor"]
 dim(tcga_data_tumot)
 tcga_data_tumot=tcga_data
 tcga_data_tumot$deceased <- ifelse(tcga_data_tumot$vital_status == "Alive", FALSE, TRUE)
 
 # create an "overall survival" variable that is equal to days_to_death
 # for dead patients, and to days_to_last_follow_up for patients who
 # are still alive
 tcga_data_tumot$overall_survival <- ifelse(tcga_data_tumot$vital_status == "Alive",
                                      tcga_data_tumot$days_to_last_follow_up,
                                      tcga_data_tumot$days_to_death)
 
 tcga_data_tumot= tcga_data_tumot [,!is.na(tcga_data_tumot$overall_survival)]
 #tcga_data_tumot= tcga_data_tumot [,colData(tcga_data_tumot)$deceased==TRUE]
 
 saveRDS(tcga_data_tumot,"E:/thesis/MSc/msc_figure/TCGA/tcga_data_tumor.rds")
 tcga_data_tumot=readRDS("E:/thesis/MSc/msc_figure/TCGA/tcga_data_tumor.rds")
 
 
 #transform=======================================
 library(DESeq2)
 #start===============

 dds <- DESeqDataSetFromMatrix(countData = assay(tcga_data_tumot),
                               colData = colData(tcga_data_tumot),
                               design = ~ 1)
 
 # Removing genes with sum total of 10 reads across all samples
 keep <- rowSums(counts(dds)) >= 10
 dds <- dds[keep,]
 
 
 # vst 
 vsd <- vst(dds, blind=FALSE)
 saveRDS(vsd,"E:/thesis/MSc/msc_figure/TCGA/tcga_data_tumor_vsd.rds")
 vsd=readRDS("E:/thesis/MSc/msc_figure/TCGA/tcga_data_tumor_vsd.rds")
 
 #======1st way========
   
   
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
 
 dge = calcNormFactors(dge,method="TMM")
 
 v = voom(dge,design,plot=TRUE)
 View(dge$sampl)
 #end=============================
 bulk.mtx=as.data.frame((v$E)) %>%
    rownames_to_column(var = 'gene_id') %>%  
    left_join(., as.data.frame(tcga_data_tumot@rowRanges)[,c("gene_id","gene_name")], by = "gene_id") 

 
 bulk.mtx$gene_id=NULL
 
 rownames(bulk.mtx)=make.unique(bulk.mtx$gene_name)
 bulk.mtx$gene_name=NULL
 bulk.mtx_t=t(bulk.mtx)
 bulk.mtx_t=rownames_to_column(as.data.frame(bulk.mtx_t),"sample_id")
clinical_data=as.data.frame(colData(tcga_data_tumot))
clinical_data=clinical_data[,c("sample_type","ajcc_pathologic_stage","gender","vital_status", "days_to_last_follow_up", "days_to_death","deceased","overall_survival")]
clinical_data=rownames_to_column(clinical_data,"sample_id")

merge_lihc=merge(clinical_data, bulk.mtx_t,by="sample_id")
dim(merge_lihc)
 View(bulk.mtx)
 
 #write.csv(merge_lihc,"E:/thesis/MSc/msc_figure/TCGA/tcga_hcc_dge.csv")
 merge_lihc=readRDS("E:/thesis/MSc/msc_figure/TCGA/tcga_hcc_vsd.rds")

 
 
 
 obj.hcc_marker=read.csv("E:/thesis/MSc/msc_figure/hcc_icca/obj.hcc_marker.csv")
 obj.hcc_marker_sig  = obj.hcc_marker  %>%
    slice_max(n = 200, order_by = (avg_log2FC))
 
 common=intersect(obj.hcc_marker_sig$X,colnames(merge_lihc))
 
 common=intersect(obj.hcc_marker_sig$X,colnames(merge_lihc))
 common=common[str_detect(common,'MT-') == FALSE]
View(obj.hcc_marker_200)

   
   
   #univariate regression==========================================================================
   
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
common=rownames(res_chol)
 covariates= common

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
   res <- t(as.data.frame(univ_results, check.names = FALSE))
   as.data.frame(res)
   write.csv(as.data.frame(res),"E:/thesis/MSc/msc_figure/TCGA/merge_lihc_univariate_Regression.csv")
   res=read.csv("E:/thesis/MSc/msc_figure/TCGA/merge_lihc_univariate_Regression.csv",row.names = 1)
   res=filter(res, p.value<0.05)
   write.csv((res),"E:/thesis/MSc/msc_figure/TCGA/coxph/merge_lihc_univariate_Regression_0.05.csv")
#end==========================================================================================================
   
   
   
#Multivariate regression======================================================================================
   library(ggplot2)
   install.packages("forestmodel")
   common_pos=paste(rownames(res), collapse="+")
   f1=as.formula(paste("Surv(overall_survival, deceased) ~ ",
                       common_pos))
   
   gene.cox <- coxph(f1,data =  (merge_lihc) )
   #suvf=survfit(f1,data =  (merge_lihc) )
   
   library(forestmodel)
   
   ggforest(gene.cox)
   
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

  
   
   forestmodel::forest_model(gene.cox,covariates = rownames(sm_cox_df),panels)
   #ggsave("E:/thesis/MSc/msc_figure/TCGA/coxph/HCC_multivariate_prognosis_forestplot.jpg", width = 8, height = 4, units = c("in"), dpi = 300)
   sm_cox=summary(gene.cox)
   sm_cox_df=as.data.frame(sm_cox$coefficients)
   
   sm_cox_df=filter(sm_cox_df, sm_cox_df$`Pr(>|z|)`<0.05)
   sm_cox_df_coeff=as.data.frame(t(select(sm_cox_df,coef)))
   rownames(sm_cox_df_coeff)=NULL
   sm_cox_df_coeff$CCL20
   
   #===============Calculating risk score from coxph significant gene==========
   
   merge_lihc$riskScore=(merge_lihc$ALDOB*sm_cox_df_coeff$ALDOB)+(merge_lihc$CHCHD10*sm_cox_df_coeff$CHCHD10)+
                        (merge_lihc$CCL20*sm_cox_df_coeff$CCL20)+(merge_lihc$LEAP2*sm_cox_df_coeff$LEAP2)+
                        (merge_lihc$SPP1*sm_cox_df_coeff$SPP1)+(merge_lihc$MIF*sm_cox_df_coeff$MIF)+
                        (merge_lihc$HMGCS2*sm_cox_df_coeff$HMGCS2)+(merge_lihc$GCHFR*sm_cox_df_coeff$GCHFR)+
                        (merge_lihc$GAMT*sm_cox_df_coeff$GAMT)+(merge_lihc$HILPDA*sm_cox_df_coeff$HILPDA)
   
   merge_lihc$Survival_years= merge_lihc$overall_survival/365
   median_value_rs = median(merge_lihc$riskScore)
   merge_lihc$RS_Group = ifelse(  merge_lihc$riskScore >= median_value_rs,  "High risk","Low risk")
      View(select(merge_lihc,riskScore,RS_Group))
      summary(fit_rs)
      table(merge_lihc$RS_Group)
      
   fit_rs = survfit(Surv(Survival_years, deceased) ~ RS_Group, data=merge_lihc)
   ggsurv=ggsurvplot(fit_rs, data=merge_lihc,surv.median.line = "none", pval=T, risk.table=F,
              xlab = "Time in years",legend.labs = 
                 c("High risk group", "Low risk group"))
   
   
   
   ggsurv$plot+ ggplot2::annotate("text", 
   x = 3, y = 0.05, # x and y coordinates of the text
   label = "High risk group (Median) - 2.75 years\nLow risk group (Median) - 6.81 years", size = 3.5)
   
   
   #ggsave("E:/thesis/MSc/msc_figure/TCGA/coxph/HCC_multivariate_prognosis.jpg", width = 5.5, height = 5, units = c("in"), dpi = 300)
#end============================================================================================================
 
 
 

 ######=========================ICCA=============================================================================
 
 
 obj.hcc_icca <- readRDS("E:/thesis/MSc/msc_figure/hcc_icca/obj.hcc_icca.rds")
 DefaultAssay(obj.hcc_icca)="RNA"
 obj.hcc_icca$orig.ident=
 #====
 obj.icc=subset(obj.hcc_icca, subset=Type=="Malignant cell")
 type=obj.icc$orig.ident
 obj.icc <- CreateSeuratObject(counts = obj.icc@assays$RNA@counts, 
                               project = "icc")
 
 type=obj.icc$orig.ident
 obj.icc=SetIdent(obj.icc, value=type)
 
 obj.icc <- NormalizeData(obj.icc)
 obj.icc  <- FindVariableFeatures(obj.icc, selection.method = "vst", nfeatures = 3000)
 all.genes <- rownames(obj.icc)
 obj.icc <- ScaleData(obj.icc, features = all.genes)
 
 
 obj.icc_marker <- FindMarkers(obj.icc, ident.1 = "iCCA",only.pos=T,logfc.threshold = 0.25)
 #====
 obj.icc=subset(obj.hcc_icca, subset=orig.ident=="iCCA")
 type=obj.icc$Type
 dim( obj.hcc_icca)
 obj.icc <- CreateSeuratObject(counts = obj.icc@assays$RNA@counts, 
                               project = "icc")
 
 
 obj.icc=SetIdent(obj.icc, value=type)
 
 obj.icc <- NormalizeData(obj.icc)
 obj.icc  <- FindVariableFeatures(obj.icc, selection.method = "vst", nfeatures = 3000)
 all.genes <- rownames(obj.icc)
 obj.icc <- ScaleData(obj.icc, features = all.genes)
 
 
 obj.icc_marker <- FindMarkers(obj.icc, ident.1 = "Malignant cell",only.pos=T,logfc.threshold = 0.25)
 obj.icc_marker$up_down=NA
 obj.icc_marker$up_down[which(obj.icc_marker$avg_log2FC>0)]='up'
 obj.icc_marker$up_down[which(obj.icc_marker$avg_log2FC<0)]='down'
 #write.csv(obj.icc_marker,"E:/thesis/MSc/msc_figure/hcc_icca/obj.icc_marker_onlyposFrom_malignant.csv")
 obj.icc_marker=read.csv("E:/thesis/MSc/msc_figure/hcc_icca/obj.icc_marker.csv")

 obj.icc_marker_200  = obj.icc_marker  %>%group_by(up_down)%>%
    slice_max(n = 200, order_by = (avg_log2FC))

 
   tcga_data_chol <- readRDS("E:/thesis/MSc/msc_figure/TCGA/tcga_data_chol.RDS") 
 assay(tcga_chol)
 View(as.data.frame(colData(tcga_data_chol)))
 dim(  hol)
 
 #======1st way========
 
 
 group = factor(colData(tcga_data_chol)$definition)
 group = relevel(group, ref="Solid Tissue Normal")
 design = model.matrix(~group)
 dge = DGEList( # creating a DGEList object
    counts=assay(tcga_data_chol),
    samples=colData(tcga_data_chol),
    genes=as.data.frame(rowData(tcga_data_chol)))
 
 # filtering
 keep = filterByExpr(dge,design) # defining which genes to keep
 dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
 rm(keep) 
 
 dge = calcNormFactors(dge,method="TMM")
 
 v = voom(dge,design,plot=TRUE)
 View(dge$samples)

 data=v$E
 #end
 
   tcga_data_col_tumor= tcga_data_chol [,colData(tcga_data_chol)$sample_type == "Primary Tumor"]
 dim(  ol_tumor)
 #  ol_tumor=  hol
   ol_tumor$deceased <- ifelse(  ol_tumor$vital_status == "Alive", FALSE, TRUE)
 
 # create an "overall survival" variable that is equal to days_to_death
 # for dead patients, and to days_to_last_follow_up for patients who
 # are still alive
   ol_tumor$overall_survival <- ifelse(  ol_tumor$vital_status == "Alive",
                                              ol_tumor$days_to_last_follow_up,
                                              ol_tumor$days_to_death)
 
   ol_tumor=   ol_tumor [,!is.na(  ol_tumor$overall_survival)]
 #  ol_tumor=   ol_tumor [,colData(  ol_tumor)$deceased==TRUE]
 
 saveRDS(  ol_tumor,"E:/thesis/MSc/msc_figure/TCGA/tcga_data_tumor.rds")
   ol_tumor=readRDS("E:/thesis/MSc/msc_figure/TCGA/tcga_data_tumor.rds")
 
 
 #transform=======================================
 library(DESeq2)
 

 #secondway===============
 
 dds <- DESeqDataSetFromMatrix(countData = assay(  ol_tumor),
                               colData = colData(  ol_tumor),
                               design = ~ 1)
 
 # Removing genes with sum total of 10 reads across all samples
 keep <- rowSums(counts(dds)) >= 10
 dds <- dds[keep,]
 
 
 # vst 
 vsd <- vst(dds, blind=FALSE)
 dim(vsd)
 saveRDS(vsd,"E:/thesis/MSc/msc_figure/TCGA/tcga_data_tumor_vsd.rds")
 vsd=readRDS("E:/thesis/MSc/msc_figure/TCGA/tcga_data_tumor_vsd.rds")
 
 #test
 
 # Get data for TP53 gene and add gene metadata information to it -------------
 gene_metadata <- as.data.frame(rowData(  ol_tumor))
 hav <- assay(vsd) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'gene_id') %>% 
    gather(key = 'case_id', value = 'counts', -gene_id) %>% 
    left_join(., gene_metadata, by = "gene_id") %>% 
    filter(gene_name == "HAVCR2")
 
 
 # get median value
 median_value <- median(hav$counts)
 
 # denote which cases have higher or lower expression than median count
 brca_tp53$strata <- ifelse(brca_tp53$counts >= median_value, "HIGH", "LOW")
 
 # Add clinical information to brca_tp53
 brca_tp53$case_id <- gsub('-01.*', '', brca_tp53$case_id)
 brca_tp53 <- merge(brca_tp53, clinical_brca, by.x = 'case_id', by.y = 'submitter_id')
 #end=============================
 bulk.mtx=as.data.frame(assay(vsd)) %>%
    rownames_to_column(var = 'gene_id') %>%  
    left_join(., as.data.frame(  ol_tumor@rowRanges)[,c("gene_id","gene_name")], by = "gene_id") 
 
 
 bulk.mtx$gene_id=NULL
 
 rownames(bulk.mtx)=make.unique(bulk.mtx$gene_name)
 bulk.mtx$gene_name=NULL
 bulk.mtx_t=t(bulk.mtx)
 bulk.mtx_t=rownames_to_column(as.data.frame(bulk.mtx_t),"sample_id")
 clinical_data=as.data.frame(colData(  ol_tumor))
 clinical_data=clinical_data[,c("sample_type","gender","vital_status","tissue_or_organ_of_origin", "days_to_last_follow_up", "days_to_death","deceased","overall_survival")]
 clinical_data=rownames_to_column(clinical_data,"sample_id")
 
 merge_lihc=merge(clinical_data, bulk.mtx_t,by="sample_id")
 dim(merge_lihc)
 View(bulk.mtx)
 #write.csv(merge_lihc,"E:/thesis/MSc/msc_figure/hcc_icca/tcga_chol.csv")
 merge_lihc=read.csv("E:/thesis/MSc/msc_figure/hcc_icca/merge_icc_dds.csv")
 
 #write.csv(merge_lihc,"E:/thesis/MSc/msc_figure/hcc_icca/merge_icc_VE.csv")
 
 obj.icc_marker=read.csv("E:/thesis/MSc/msc_figure/hcc_icca/obj.icc_marker_onlyposFrom_malignant.csv")
 dim( obj.icc_marker)
 obj.icc_marker= obj.icc_marker[obj.icc_marker$up_down=='up',]
 obj.icc_marker_sig =filter(obj.icc_marker, avg_log2FC<=0.5)
 obj.icc_marker_sig  = obj.icc_marker %>%
    slice_max(n = 1200, order_by = (avg_log2FC))
 
 common=intersect(obj.icc_marker_sig$X,colnames(merge_lihc))
 
 common=common[str_detect(common,'MT-') == FALSE]
 common=common[str_detect(common,'MIR4435-2HG|AS1|HLA|PHGR1|ERV3') == FALSE]
 View(common)
 
 MGST1
 coxph(Surv(overall_survival, deceased)~ MGST1,data = merge_lihc)
 
 #univariate regression==========================================================================
 
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
 
 covariates= common
 
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
                       
 write.csv(as.data.frame(res_chol),"E:/thesis/MSc/msc_figure/TCGA/merge_icc_univariate_Regression.csv")
 res_chol=read.csv("E:/thesis/MSc/msc_figure/TCGA/merge_icc_univariate_Regression.csv",row.names = 1)
 res_chol=filter(as.data.frame(res_chol), p.value<=0.05)
 
 filter(obj.icc_marker, X %in% rownames(res_chol))
 rownames(res_chol)
 write.csv((res_chol),"E:/thesis/MSc/msc_figure/TCGA/coxph/merge_icc_univariate_Regression_0.05.csv")
 #end==========================================================================================================
 
 
 
 #Multivariate regression======================================================================================
 library(ggplot2)
 install.packages("forestmodel")
 common_pos=paste(rownames(res_chol), collapse="+")
 f1=as.formula(paste("Surv(overall_survival, deceased) ~ ",
                     common_pos))
 
 gene.cox <- coxph(f1,data =  (merge_lihc) )
 #suvf=survfit(f1,data =  (merge_lihc) )
 sm_cox=summary(gene.cox)
 sm_cox
 sm_cox_df=as.data.frame(sm_cox$coefficients)
 
 sm_cox_df=filter(sm_cox_df, sm_cox_df$`Pr(>|z|)`<=0.05)
 sm_cox_df=filter(sm_cox_df, sm_cox_df$`Pr(>|z|)`<=0.05)
 library(forestmodel)
 
 ggforest(gene.cox)
 
 panels <- list(
    list(width = 0.03),
    list(width = 0.09,fontsize=20, display = ~variable, fontface = "bold", heading = "Gene"),
    
    
    
    list(width = 0.03, item = "vline", hjust = 0.5),
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
 
 
 
 forestmodel::forest_model(gene.cox,covariates = rownames(sm_cox_df),panels)
 ggsave("E:/thesis/MSc/msc_figure/TCGA/coxph/icca_multivariate_prognosis.jpg", width = 5, height = 4, units = c("in"), dpi = 300)

 sm_cox_df_coeff=as.data.frame(t(select(sm_cox_df,coef)))
 rownames(sm_cox_df_coeff)=NULL
 sm_cox_df_coeff$CCL20
 
 #===========AIC test=================
 install.packages("AICcmodavg")
 library(AICcmodavg)
 
 
 interaction.mod <- lm(overall_survival ~ HMGA1+KLK6+CDH1+TMEM45A+CFH+DCBLD2+EGFR+STEAP1+SDC1+KYNU+PYGB+PLOD2+ST14+SH3BGRL2, data = merge_lihc)
 
 #===============Calculating risk score from coxph significant gene==========
 
 merge_lihc$riskScore=(merge_lihc$PLCG2*sm_cox_df_coeff$PLCG2)+(merge_lihc$SH3BGRL2*sm_cox_df_coeff$SH3BGRL2)+
    (merge_lihc$C20orf27*sm_cox_df_coeff$C20orf27)+(merge_lihc$HMGA2*sm_cox_df_coeff$HMGA2)
 
 merge_lihc$Survival_years= merge_lihc$overall_survival/365
 median_value_rs = median(merge_lihc$riskScore)
 merge_lihc$RS_Group = ifelse(  merge_lihc$riskScore >= median_value_rs,  "High risk","Low risk")
 View(select(merge_lihc,riskScore,RS_Group))
 summary(fit_rs)
 table(merge_lihc$RS_Group)
 
 fit_rs = survfit(Surv(Survival_years, deceased) ~ RS_Group, data=merge_lihc)
 ggsurv=ggsurvplot(fit_rs, data=merge_lihc,surv.median.line = "none", pval=T, risk.table=F,
                   xlab = "Time in years",legend.labs = 
                      c("High risk group", "Low risk group"))
 
 
 
 ggsurv$plot+ ggplot2::annotate("text", 
                                x = 1.6, y = 0.07, # x and y coordinates of the text
                                label = "High risk group (Median) - 1.22 years\nLow risk group (Median) - 5.31 years", size = 3.5)
 
 
 
 #========================Go pathway analyis===============
 
 SYMBOL2EG <-
    eval(parse(text = sprintf(
       'org.%s.egSYMBOL2EG', 'Hs'
    )))
 logFC_score <- obj.icc_marker$avg_log2FC
 logFC_score
 de_genes<- (obj.icc_marker$X)
 length(x = de_genes)
 names(logFC_score) <- de_genes
 logFC_score
 
 genes <- intersect(de_genes, mappedkeys(SYMBOL2EG)) #  Get the gene symbol that are mapped to an entrez gene identifiers
 #length(x = genes)
 
 
 logFC_score<- logFC_score[genes] # access logFC score 
 #length(x = logFC_score) 
 gene_entrez <-genes %>% SYMBOL2EG[.] %>% as.list %>% map( ~ .[1]) %>% purrr::simplify()
 #length(x = gene_entrez)
 
 names(logFC_score) <- gene_entrez
 logFC_score # input for rank in fgsea
 
 
 
 
 
 gse <- gseGO(geneList=sort(logFC_score, decreasing = T), 
              ont ="BP", 
              scoreType = "pos",
              nPermSimple=5000,
              
              OrgDb = org.Hs.eg.db
 )
 
 View(gse@result)
 remove.packages("clusterProfiler")
 
# saveRDS(gse,"E:/thesis/MSc/msc_figure/hcc_icca/gse_hcc_malignant.rds")
 gse=readRDS("E:/thesis/MSc/msc_figure/hcc_icca/gse_hcc_malignant.rds")
 gse1=filter(gse, !str_detect(Description,"positive|tube|blood|protein|plasma"))
 gse1=(gse1%>% arrange(desc(enrichmentScore)))
 gse1_20=gse1[1:20]
 
 dotplot(gse1, showCategory=20) 
 barplot(gse1, showCategory=20) 
 barplot((gse1), showCategory=20) 
 dotplot(gse1, x="enrichmentScore", showCategory=20)
 
 
 
 
 
 
 gse1_20$GeneRatio=(round(-log10(gse1_20$enrichmentScore),2))
 
 gse1_20$GeneRatio= gse1_20$enrichmentScore
 gse1_20=gse1_20%>%arrange(GeneRatio,decreasing=T)
 typeof(gse1_20$GeneRatio)
 ggplot(gse1_20, 
        aes(x=reorder(Description,enrichmentScore), 
            y= enrichmentScore,
            fill="Upregulated")) +  labs(title = paste("GO"))+
    
    geom_bar(stat = "identity",width = 0.6, show.legend = FALSE)+coord_flip()+
    theme(panel.background = element_blank(),
          axis.line.x.top = element_line(size = 1),
          axis.line.x.bottom = element_line(size = 0.3),
          axis.line.y.left = element_line(linetype="dotted"),
          axis.text=element_text(color="black",size=10),
          axis.title.y =element_blank(),
          plot.title = element_text(color="black",size=10,face = "bold",hjust = 0.5) )
 
 
 #ggsave("E:/thesis/MSc/msc_figure/hcc_icca/icc_malignantcell_go_enrichment.jpg", width = 6.5, height = 4, units = c("in"), dpi = 300)
 
 
 
 gsekegg <- gseKEGG(geneList     = sort(logFC_score, decreasing = T),
                    organism     = 'hsa',
                    scoreType = "pos",
                    pvalueCutoff = 1,
                    verbose      = FALSE)
 View(gsekegg1@result)
 
# saveRDS(gsekegg,"E:/thesis/MSc/msc_figure/hcc_icca/gsekegg_icc_malignant.rds")
 #gsekegg=readRDS("E:/thesis/MSc/msc_figure/hcc_icca/gsekegg_hcc_malignant.rds")
 gsekegg1=filter(gsekegg, !str_detect(Description,"Coronavirus|lung|Amoebiasis|cardiomyopathy"))
 gsekegg1=(gsekegg1%>% arrange(desc(enrichmentScore)))
 gsekegg1_20=gsekegg1[1:20]
 
 library(ggplot)
 
 dotplot(gsekegg, showCategory=20)  + scale_color_gradient("Size of\nthe pathway" , low="blue", high="red") 
 #write.csv( gsekegg1_20,"E:/thesis/MSc/msc_figure/hcc_icca/gsekegg_hcc20.csv")
 
 #ggsave("E:/thesis/MSc/msc_figure/hcc_icca/icc_malignantcell_KEGG_enrichment.jpg", width = 6, height = 4, units = c("in"), dpi = 300)
 gsekegg1_20$GeneRatio=(round(-log10(gsekegg1_20$enrichmentScore),2))
 gsekegg1_20=gsekegg1_20%>%arrange(GeneRatio,decreasing=T)
 gsekegg1_20$GeneRatio= gsekegg1_20$enrichmentScore
 typeof(gsekegg1_20$GeneRatio)
 ggplot(gsekegg1_20, 
        aes(x=reorder(Description,GeneRatio), 
            y=GeneRatio,
            fill="Upregulated")) +  labs(title = paste("KEGG"))+
    geom_bar(stat = "identity",width = 0.6, show.legend = FALSE)+coord_flip()+
    theme(panel.background = element_blank(),
          axis.line.x.top = element_line(size = 1),
          axis.line.x.bottom = element_line(size = 0.3),
          axis.line.y.left = element_line(linetype="dotted"),
          axis.text=element_text(color="black",size=10),
          axis.title.y =element_blank(),
          plot.title = element_text(color="black",size=10,face = "bold",hjust = 0.5) )
 
 
 
 
 
 
 #====================gene ics analsis tcga chol
 merge_lihc=read.csv("E:/thesis/MSc/msc_figure/hcc_icca/merge_icc_dds.csv")
 img=read.csv("E:/thesis/MSc/msc_figure/TCGA/InnateDB_genes.csv")
 img_gene=img$name
 tcga_chol=read.csv("E:/thesis/MSc/msc_figure/hcc_icca/tcga_chol.csv")
 
 tcga_chol_icc=tcga_chol[which(tcga_chol$tissue_or_organ_of_origin=="Intrahepatic bile duct"),]
 tcga_chol_icc=tcga_chol[which(tcga_chol_icc$sample_type=="Primary Tumor"),]
 dim(tcga_chol_icc)
 
 
 ics=ics=c("TIGIT", "LY9", "PDCD1", 'LAG3', "CTLA4", "CD276", "NT5E", "PDCD1LG2", "CD274", "IDO1", "VSIR", "HAVCR2", "ENTPD", "CD274", "CTLA4", "HAVCR2", "LAG3", "TIGIT",
           "PDCD1", "PDCD1LG2","SIGLEC15","VSIR","ICOS","CD80","TNFRSF18","TNFRSF4","TNFRSF9",
           "ISYNA1","CSF1","CD47","TLR8","ARG1","TREM2BMI1","TNFRSF1B","BTK","LILRB1","PLK1","CD52","CCR5",
           "TLR3","TLR7","CREBRF","TLR9","TLR4","LILR2","CXCR4","CCR2","NOX4","CD96","CD38","AXL","NR2F6",
           "ICOS", "VISTA", "CD112R", "BTLA", "GITR", "NKG2A")
 
 ics=unique(ics)
 
 
 length(obj.hcc_marker$X)
 
 merge_lihc=tcga_chol_icc
 common=intersect(img_gene, obj.hcc_marker$X)
 common=intersect(common,colnames(merge_lihc))
 
 covariates= common
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
 Cluste1MCC_top38_default_node <- read_csv("E:/thesis/MSc/msc_figure/cytoscape/Cluste1MCC_top38 default node.csv")
 cytoscape_gene=Cluste1MCC_top38_default_node$`query term`
 res_chol=filter(as.data.frame(res_chol), rownames(res_chol)%in%cytoscape_gene)
 res_chol=filter(as.data.frame(res_chol), p.value<0.05)
 res_chol=rownames_to_column(res_chol,"gene")
 filter(obj.icc_marker, X %in% rownames(res_chol))
 rownames(res_chol)
 #write.csv(res_chol,"E:/thesis/MSc/msc_figure/cytoscape/cytoscapeGeneImmune.csv")
 library(glmnet)
x=merge_lihc[, (res_chol$gene)]
y=merge_lihc[,"overall_survival"]
rownames(merge_lihc)=merge_lihc$X
cv_model <- cv.glmnet(x, y, alpha = 1)

 #end==========================================================================================================
 
 
 
 #Multivariate regression======================================================================================
 library(ggplot2)
 install.packages("forestmodel")
 common_pos=paste(rownames(res_chol), collapse="+")
 f1=as.formula(paste("Surv(overall_survival, deceased) ~ ",
                     common_pos))
 
 gene.cox <- coxph(f1,data =  (merge_lihc) )
 #suvf=survfit(f1,data =  (merge_lihc) )
 sm_cox=summary(gene.cox)
 sm_cox
 sm_cox_df=as.data.frame(sm_cox$coefficients)
 
 sm_cox_df=filter(sm_cox_df, sm_cox_df$`Pr(>|z|)`<=0.05)
 sm_cox_df=filter(sm_cox_df, sm_cox_df$`Pr(>|z|)`<=0.05)
 library(forestmodel)
 
 ggforest(gene.cox)
 
 panels <- list(
    list(width = 0.03),
    list(width = 0.09,fontsize=20, display = ~variable, fontface = "bold", heading = "Gene"),
    
    
    
    list(width = 0.03, item = "vline", hjust = 0.5),
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
 
 
 
 forestmodel::forest_model(gene.cox,covariates = rownames(sm_cox_df),panels)
 
 #===========================================================
 # get median value
 HAVCR2_median_value <- median(tcga_chol_icc$HAVCR2)
 
 # denote which cases have higher or lower expression than median count
 tcga_chol_icc$HAVCR2_HL <- ifelse(tcga_chol_icc$HAVCR2 >= HAVCR2_median_value, "HIGH", "LOW")
 
 coxph(Surv(overall_survival, deceased)~IGF2, tcga_chol) 
 
 # fitting survival curve -----------
 fit <- survfit(Surv(overall_survival, deceased) ~ kHAVCR2_HL, data = tcga_chol_icc)
 fit
 ggsurvplot(fit,
            data = tcga_chol_icc,
            pval = T,
            risk.table = T)
 
 
 fit2 <- survdiff(Surv(overall_survival, deceased) ~ HAVCR2_HL, data = tcga_chol)
 