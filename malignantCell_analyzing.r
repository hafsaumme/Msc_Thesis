# Retrieving malignant cells from HCC samples
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

# malignant cells markers

obj.hcc_marker <- FindMarkers(obj.hcc, ident.1 = "Malignant cell",
                              only.pos=T,test.use = "wilcox", 
                              logfc.threshold = 0.25)



#========================= GO and KEGG pathway analysis===============================

# Preparing marker genes 

SYMBOL2EG <-
  eval(parse(text = sprintf(
    'org.%s.egSYMBOL2EG', 'Hs'
  )))

logFC_score <- obj.hcc_marker$avg_log2FC
logFC_score
de_genes<- rownames(obj.hcc_marker)
length(x = de_genes)
names(logFC_score) <- de_genes
logFC_score

genes <- intersect(de_genes, mappedkeys(SYMBOL2EG)) #  Get the gene symbol that are mapped to an entrez gene identifiers

ogFC_score<- logFC_score[genes] # access logFC score

gene_entrez <-genes %>% SYMBOL2EG[.] %>% as.list %>% map( ~ .[1]) %>% purrr::simplify()

names(logFC_score) <- gene_entrez

#=== GO pathways==

gse_go<- gseGO(geneList=sort(logFC_score, decreasing = T), 
             ont ="BP", 
             pvalueCutoff = 0.05,
             scoreType = "pos",
             OrgDb = org.Hs.eg.db
)

gse_go=gse_go@result
gse_go=(gse_go%>% arrange(desc(enrichmentScore)))
gse_go=gse_go[1:20,]
gse_go$GeneRatio= gse_go$enrichmentScore


ggplot(gse_go, 
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


#========KEGG pathways==========

gse_kegg <- gseKEGG(geneList  = sort(logFC_score, decreasing = T),
                   organism   = 'hsa',
                   scoreType = "pos",
                   minGSSize = 1,
                   maxGSSize = Inf,
                   pvalueCutoff = 5,
                   verbose      = FALSE)

gsekegg=gse_kegg@result
gsekegg=filter(gsekegg, pvalue<0.02)
gsekegg=(gsekegg%>% arrange(desc(enrichmentScore)))
gsekegg=gsekegg[1:12,]
gsekegg$GeneRatio= gsekegg$enrichmentScore


ggplot(gsekegg, 
       aes(x=reorder(Description,GeneRatio), 
           y= GeneRatio,
           fill = pvalue
       )) +
  geom_col(width = 0.8)+coord_flip()+ 
  scale_fill_gradient(low = "red", high = "blue")+theme_bw()+
  scale_x_discrete(labels = function(x) str_wrap(x, width =70))+
  theme(axis.text=element_text(color="black",size=12),
        panel.border = element_rect(color="black",fill=NA,linewidth =0.5),
  )+labs(title = "",
       x="",
       y = 'GeneRatio')







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

obj.hcc_marker_sig  = obj.hcc_marker  %>%slice_max(n = 200, order_by = (avg_log2FC))

common=intersect(obj.hcc_marker_sig$X,colnames(merge_lihc))

common=intersect(obj.hcc_marker_sig$X,colnames(merge_lihc))

common=common[str_detect(common,'MT-') == FALSE]

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

# selecting those genes that have pvalue<0.05
res=filter(res, p.value<0.05)


#Multivariate regression======================================================================================

common_pos=paste(rownames(res), collapse="+")
f1=as.formula(paste("Surv(overall_survival, deceased) ~ ",
                    common_pos))

gene.cox <- coxph(f1,data =  (merge_lihc) )
sm_cox=summary(gene.cox)
sm_cox
sm_cox_df=as.data.frame(sm_cox$coefficients)
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




#==============Survival analysis==========================================================================



tumor_data=merge_lihc[(merge_lihc$sample_type == "Primary Tumor"),]

# kaplan-Meier Curves for 13 hub genes from regression analyis

gene="CPB2"

median_value <- median(tumor_data[, gene])

tumor_data[,gene] <- ifelse(tumor_data[,gene] >= median_value, "HIGH", "LOW")

tumor_data$Survival_years= tumor_data$overall_survival/365

fit <-survfit(as.formula(paste('Surv(Survival_years, deceased)~', gene)),data=tumor_data)

ggsurvplot(fit, data=tumor_data,surv.median.line = "none", pval=T, risk.table=F,
           xlab = "Time in years",legend.title="",legend = c(0.8, 0.9),legend.labs = 
             c("High expression", "Low expression"))$plot+ggtitle(gene)+
  theme(plot.title = element_text(hjust=0.5, face = "bold",size=13))


# Boxplot of  genes between normal and tumor samples

ggboxplot(merge_lihc, 
          x = "sample_type", 
          y = "CPB2", 
          fill="sample_type",
          outlier.shape = NA,
          bxp.errorbar = T,
          bxp.errorbar.width = 0.2,
          xlab=FALSE,ylab="log2(TPM + 1)",title = "CPB2")+ 
  theme(legend.position = "none",plot.title = element_text(hjust=0.5, face = "bold",size=14))+
  geom_jitter(position = position_jitter(height = 4, width = .08),size=0.6)+
  geom_pwc(method = "t_test",label="p.format",bracket.nudge.y = 0.4, size=0.5,tip.length = 0.04)




#correlation analysis===========================


immunecell=(read.csv("cibersortxResult.csv",row.names = 1))
merge_lihc_hub=merge_lihc[,c("sample_id","C3","GC","CPB2","FGB","HRG","KNG1","SERPINC1")]

library(stringr)

immunecell=rownames_to_column(immunecell,"sample_id")
immunecell$sample_id=gsub("\\.", "-", immunecell$sample_id)
merge_immunecell=merge(immunecell,merge_lihc_hub,by="sample_id")
merge_immunecell=merge_immunecell%>%select(-c("P.value","Correlation","RMSE","sample_id"))


result = cor(merge_immunecell,method = "spearman")
result =result[,1:22]
result=result[23:29,]
result=round(result,2)



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
