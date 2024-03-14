setwd("E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets")

getwd()
unloadNamespace("Matrix")
install.packages('Matrix')
install.packages('Matrix', repos='https://CRAN.R-project.org/package=Matrix', dependencies=TRUE)
remove.packages(grep("spatstat", installed.packages(), value = T))
install.packages(c("spatstat.utils","spatstat.core", "spatstat.data"), 
                 repos = "https://spatstat.r-universe.dev")
package_version("seuratObject")
remove.packages("SeuratObject")
devtools::install_version("spatstat", version = "1.64-1")
install.packages('spatstat.data')
install.packages('spatstat.core')
remove.packages("Seurat")
remove.packages("Matrix")
install.packages('HGNChelper')
install.packages("SeuratObject")
remotes::install_version("Seurat", "3")
Biobase::package.version("Matrix")
ibrary(dplyr)
options("install.lock"=FALSE)
#library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(tibble)
library(dplyr)
library(stringr)
library(glmGamPoi)
#install.packages("Matrix", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
#p_unlock("C:/Program Files/R/R-4.3.1/library")
#p_unlock(lib.loc = p_path("Matrix"))
Sys.getenv()
library(pacman)
install.packages("pacman")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("glmGamPoi")
obj.plasma_cirrhotic_healthy= readRDS("E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/obj.plasma_cirrhotic_healthy.rds")
dim(obj.plasma_cirrhotic_healthy)# 23039   913
obj.plasma_cirrhotic_healthy$sample=obj.plasma_cirrhotic_healthy$Liver
obj.plasma_cirrhotic_healthy <- SCTransform(obj.plasma_cirrhotic_healthy, method = "glmGamPoi", 
                                 vars.to.regress = "percent.mt",variable.features.n = 3000,
                                 verbose = FALSE)

ob.bcell.hcc= (readRDS("E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/obj.bcell.HCC.rds"))
ob.bcell.hcc$sample=ob.bcell.hcc$liver
DefaultAssay(ob.bcell.hcc) <- "RNA"
ob.bcell.hcc <- SCTransform(ob.bcell.hcc, method = "glmGamPoi", 
                            vars.to.regress = "percent.mt",variable.features.n = 3000,
                            verbose = FALSE)
dim(ob.bcell.hcc)#17814   605


ob.bcell.icca =(readRDS("materials/hafsa/uh_msc/gsc136103/obj.bcell.icca.rds"))
ob.bcell.icca$sample=ob.bcell.icca$liver
DefaultAssay(ob.bcell.icca) <- "RNA"
ob.bcell.icca <- SCTransform(ob.bcell.icca, method = "glmGamPoi", 
                            vars.to.regress = "percent.mt",variable.features.n = 3000,
                            verbose = FALSE)
dim(ob.bcell.icca) #18009   388

ob.bcell.cirrhotic_healthy =(readRDS("materials/hafsa/uh_msc/gsc136103/obj.bcell_cirrhotic_healthy.rds"))
ob.bcell.cirrhotic_healthy$sample=ob.bcell.cirrhotic_healthy$Liver

DefaultAssay(ob.bcell.cirrhotic_healthy) <- "RNA"
dim(ob.bcell.cirrhotic_healthy) #23039  1929

obj.bcell.cirrhotic<- subset(ob.bcell.cirrhotic_healthy, subset = orig.ident == 'cirrhotic')
dim(obj.bcell.cirrhotic) #21085  1157
obj.bcell.cirrhotic <- SCTransform(obj.bcell.cirrhotic, method = "glmGamPoi", 
                            vars.to.regress = "percent.mt",variable.features.n = 3000,
                            verbose = FALSE)

obj.bcell.healthy<- subset(ob.bcell.cirrhotic_healthy, subset = orig.ident == 'healthy')
dim(obj.bcell.healthy) #21085   772
obj.bcell.healthy <- SCTransform(obj.bcell.healthy, method = "glmGamPoi", 
                                   vars.to.regress = "percent.mt",variable.features.n = 3000,
                                   verbose = FALSE)


#....................INTEGRATION on hcc , icca, healthy and cirrhotic datasets...............

list= c(obj.plasma_cirrhotic_healthy, obj.bcell.healthy,obj.bcell.cirrhotic,
        ob.bcell.icca,ob.bcell.hcc)


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 4000)

list <- PrepSCTIntegration(object.list = list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features)
obj.all.four <- IntegrateData(anchorset = anchors,normalization.method = "SCT")

DefaultAssay(obj.all.four ) <- "RNA"
dim(obj.all.four@assays$integrated@var.features)#4000 3835
View(obj.all.four@meta.data)



obj.all.four@meta.data$Liver=NULL
obj.all.four@meta.data$liver=NULL
obj.all.four@meta.data$integrated_snn_res.0.5=NULL
obj.all.four@meta.data$seurat_clusters=NULL
obj.all.four@meta.data$newCellName=NULL

obj.all.four@meta.data$type[which(str_detect(obj.all.four@meta.data$orig.ident,"iCCA"))] <- 'tumor'
obj.all.four@meta.data$type[which(str_detect(obj.all.four@meta.data$orig.ident,"healthy"))] <- 'Healthy'
obj.all.four@meta.data$customclassif[which(str_detect(obj.all.four@meta.data$Type,"B cell"))] <- 'Bcell'
obj.all.four@meta.data$Type=NULL
#View(obj.all.four@meta.data)

#Dimensionality Reduction



obj.all.four=RunPCA(obj.all.four, features = VariableFeatures(object = obj.all.four))

ElbowPlot(obj.all.four, ndims = 50)

pct <- obj.all.four[["pca"]]@stdev / sum(obj.all.four[["pca"]]@stdev) * 100

cumu <- cumsum(pct)

co1 <- which(cumu > 90 & pct < 5)[1]

co1
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2
pcs <- min(co1, co2)

pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

library(ggplot2)
# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()



obj.all.four=RunUMAP(obj.all.four, reduction = "pca", dims = 1:pcs)
dim(obj.all.four)
DimPlot(obj.all.four, reduction = "umap")
obj.all.four@meta.data=NULL
#obj.all.four=ob.bcell.hcc
obj.all.four <- FindNeighbors(obj.all.four, reduction = "pca", dims = 1:15)
obj.all.four= FindClusters(obj.all.four,resolution = c(0.1))
obj.all.four1=(readRDS("obj.all.four.bcell1.rds"))
obj.all.four1 <- FindNeighbors(obj.all.four1, reduction = "pca", dims = 1:15)
obj.all.four1= FindClusters(obj.all.four1,resolution = c(0.4))
obj.all.four1=RunUMAP(obj.all.four1, reduction = "pca", dims = 1:15)
DimPlot(obj.all.four1, reduction = "umap")
#install.packages("clustree")
#install.packages("igraph")
#install.packages("igraph")
library(igraph)
library(clustree)

clustree::clustree(obj.all.four1)
obj.all.four=obj.all.four1


obj.all.four@assays$RNA@data=NormalizeData(obj.all.four@assays$RNA@counts, normalization.method = "LogNormalize", scale.factor = 10000)
obj.all.four.genes <- rownames(obj.all.four)
obj.all.four@assays$RNA@scale.data= ScaleData(obj.all.four@assays$RNA@data, features = rownames(obj.all.four))
DefaultAssay(obj.all.four ) <- "integrated"

obj.all.four_marker=FindAllMarkers(obj.all.four, test.use = "wilcox", logfc.threshold = 0.25)
table(obj.all.four$orig.ident)

obj.all.four_marker = obj.all.four_marker %>%
  group_by(cluster) %>%
  slice_max(n = 500, order_by = avg_log2FC)
#View(obj.all.four_marker)

#==============sctype

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "E:/thesis/MSc/msc data/ScTypeDB_full.xlsx"; tissue = "Liver_Bcell"
#View(db_)
#db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"; 









gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = obj.all.four@assays$integrated@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
#View(es.max)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(obj.all.four@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(obj.all.four@meta.data[obj.all.four@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(obj.all.four@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
#(View(sctype_scores))
# set low-confident (low ScType score) clusters to "unknown"
#sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
#print(sctype_scores[,1:3])



obj.all.four@meta.data$subtype = NA
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  obj.all.four@meta.data$subtype[obj.all.four@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

#View(combined@meta.data)
DimPlot(obj.all.four, reduction = "umap", label = F, repel = T)+
  theme(axis.title = element_text(size=11))+ggtitle("")+labs(color="")
#saveRDS(HSMM2, "E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/trajectory/HSMM2_bcell.rds")
ggsave("E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/bcell_type_seurat.jpg", dpi=300)
View(obj.all.four@meta.data)
#saveRDS(obj.all.four,"materials/hafsa/uh_msc/bcell_from_all_datasets/obj.all.four.bcell.rds")
obj.all.four=(readRDS("obj.all.four.bcell.rds"))
obj.all.four1=subset(obj.all.four, subset=(condition=="Tumor"|condition=="Healthy"))
obj.all.four1=subset(obj.all.four, subset=(orig.ident=="HCC"| orig.ident=="healthy"))
obj.all.four1=SetIdent(obj.all.four1,value=obj.all.four1$condition)
obj.all.four1=subset(obj.all.four1, subset=customclassif=="Naïve B cells")
obj.all.four1_count=obj.all.four1@assays$RNA@counts
obj.all.four1_count=obj.all.four1_count[obj.all.four_marker$gene,] 
colnames(obj.all.four1_count)=obj.all.four1$customclassif
colnames(obj.all.four1_count)=sub("\\.\\d+$","",colnames(obj.all.four1_count))
obj.all.four1_count=obj.all.four1_count[,order(names(obj.all.four1_count))]
obj.all.four1_count=as.data.frame(obj.all.four1_count)
colnames(obj.all.four1_count)="Naive.Bcell"
dim(obj.all.four1_count)
obj.all.four1_meta=obj.all.four1@meta.data%>% group_by(customclassif, condition) %>% 
  summarise(frequency = n()) %>%
  mutate(C = sum(frequency)) %>%
  mutate(percent = frequency/C*100)
#unloadNamespace("clusterprofile")
ggplot(obj.all.four1_meta, aes(x = "", y = "percent", fill = "condition")) +
  geom_col() +
  coord_polar(theta = "y") +
  geom_text(aes(label = paste(percent, "%")),
            position = position_stack(vjust = 0.5),
            size = 8) +
  theme_void(base_size = 20) +
  scale_fill_brewer(name = NULL, palette = "Pastel2")
display.brewer.all()
obj.all.four1_meta$percent=round(obj.all.four1_meta$percent,2)
ggplot(obj.all.four1_meta, aes(x="", y=percent, fill=condition)) +labs(title = "Naive B cells")+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(label = paste(percent, "%")),
            position = position_stack(vjust = 0.5),
            
            size = 5) + scale_fill_brewer(palette="Set2")+
  
  theme_void()+theme(plot.title=element_text(hjust=0.5,vjust=-4,face="bold"))

DotPlot(object=obj.all.four,features = c("CD79A","IGHM","CD79B","CD37","CD52","MS4A1","CD27","CD38","CD40","JCHAIN",
                                         "IGHA1","IGHG4","IGHA2","IGHG3","IGHG2"))+
                                          theme(axis.text.x = element_text(angle=60,hjust=1))
ggsave("E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/hcc_icc_type.jpg",
       width = 5, height = 3.4, units = c("in"), dpi = 300)
obj.all.four@meta.data$customclassif[which(rownames(obj.all.four@meta.data) ==colnames(obj.plasma_cirrhotic_healthy))]="Plasma cells"
obj.all.four@meta.data$type=NA
obj.all.four@meta.data$type[which(str_detect(rownames(obj.all.four@meta.data), rownames(obj.plasma_cirrhotic_healthy@meta.data))),] = "Plasma cell"

obj.all.four@meta.data$type[which(is.na(obj.all.four@meta.data$type))]="B cell"
str_detect("YWnjjj", "[YW]")
View(obj.all.four@meta.data)
# Rename all identities
obj.all.four <- RenameIdents(object = obj.all.four, 
                         "0" = "Naive Bcells",
                         "1" = "Naive Bcells",
                         "2" = "Memory Bcells",
                         "3" = "Plasma cells",
                         "4" = "Memory Bcells",
                         "5" = "Plasma cells",
                         "6" = "Memory Bcells"
                         )

table(obj.all.four$customclassif)
DimPlot(obj.all.four, reduction = "umap", group.by ="customclassif",label=T)


DefaultAssay(obj.all.four) <- "RNA"
obj.all.four=SetIdent(obj.all.four,value=obj.all.four$subtype)
obj.all.four_marker=FindAllMarkers(obj.all.four, only.pos = TRUE, test.use = "wilcox", logfc.threshold = 0.25)
table(obj.all.four$customclassif)
#View(obj.all.four_marker)

#write.csv(obj.all.four1_count,"E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/signaturenaivebcell3.csv")
obj.all.four_marker=((read.csv("obj.all.four_marker.csv")))
obj.all.four_marker=obj.all.four_marker[which(obj.all.four_marker$cluster=="Naïve Bcells"),]
obj.all.four=SetIdent(obj.all.four,value = obj.all.four$orig.ident)
View(obj.all.four_marker)
VlnPlot(object = obj.all.four1, 
        features = c("MS4A1"))
obj.all.four_marker = obj.all.four_marker %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

obj.all.four_marker=marker_orig.ident
##end

#obj.all.four=SCTransform(obj.all.four, method = "glmGamPoi", 
           # vars.to.regress = "percent.mt",variable.features.n = 3000,
            #verbose = FALSE)
obj.all.four@assays$RNA@data=NormalizeData(obj.all.four@assays$RNA@counts, normalization.method = "LogNormalize", scale.factor = 10000)
obj.all.four <- FindVariableFeatures(obj.all.four, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(obj.all.four)
obj.all.four@assays$RNA@scale.data= ScaleData(obj.all.four@assays$RNA@data, features = all.genes)


dim(obj.all.four@assays$RNA)
obj.all.four@meta.data$condition[obj.all.four@meta.data$orig.ident == "healthy"] = "Healthy"
obj.all.four@meta.data$condition[obj.all.four@meta.data$orig.ident == "iCCA"] = "Tumor"
o
obj.all.four<- SetIdent(obj.all.four, value = obj.all.four@meta.data$orig.ident)

marker_orig.ident=FindAllMarkers(obj.all.four, only.pos = TRUE, test.use = "wilcox", logfc.threshold = 0.25)
marker_orig.ident=(read.csv("orig.ident_marker.csv"))
marker_orig.ident = marker_orig.ident %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)
View(obj.all.four@meta.data)
#trajectory analysis

hcc.markers <- FindMarkers(obj.all.four, ident.1 = "HCC", only.pos = TRUE, test.use = "wilcox", logfc.threshold = 0.25)
hcc.markers= hcc.markers %>%
  slice_max(n = 20, order_by = avg_log2FC)
head(cluster5.markers, n = 5)



library(monocle)
#devtools::install_github("cysouw/qlcMatrix")
#BiocManager::install("monocle")

data <- as(as.matrix(obj.all.four@assays$RNA@scale.data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data =obj.all.four@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
#fData

fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
HSMM2 <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        #lowerDetectionLimit = 0.5,
                        expressionFamily =uninormal())


# Select genes (High varible gene/all gene) and Run ordering algorithm for pseudotime
#ordering_genes<-(obj.all.four)@assays$integrated@var.features
ordering_genes<-(obj.all.four_marker$X)
#vargene=(c("CD27", "CD79A","IGHA1", "IGHA2", "CD27","CD38","TCL1A", "CD19","MS4A1") )

#vargene=ordering_genes
#ordering_genes<-c("CD27", "CD79A","IGHA1", "IGHA2", "CD27","CD38","TCL1A", "CD19","MS4A1")
length(ordering_genes)
#ordering_genes<- rownames(x = combined.data[["RNA"]]) #for all gene


HSMM2 <- setOrderingFilter(HSMM2, ordering_genes)

## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
HSMM2 <- reduceDimension(HSMM2, norm_method="none",
                         reduction_method="DDRTree",
                         max_components=2,
                         scaling=TRUE,
                         verbose=TRUE,
                         pseudo_expr=0)
#View data
pData(HSMM2)$celltype=obj.all.four@meta.data$customclassif
#pData(HSMM2)$Source=pData(HSMM2)$orig.ident# Attribute of cell
pData(HSMM2)
fData(HSMM2) # gene attribute

HSMM2 <- orderCells(HSMM2,root_state=3)
HSMM2


# .....trajectory plot
HSMM2=readRDS("E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/trajectory/HSMM2_bcell.rds")

HSMM_myo <- setOrderingFilter(HSMM2, ordering_genes)
plot_ordering_genes(HSMM_myo)
# Three option for 'color_by'= "seurat_clusters"/"Pseudotime"/"State"
plot_cell_trajectory(HSMM2, 
                     color_by = "Source",
                     theta = -15,
                     show_branch_points = T,
                     show_tree = TRUE,
                     cell_size =1)+
  
  facet_wrap(~Source, nrow = 1)

ggsave("E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/ms4a1_inscrna.jpg",
       width = 4, height = 4, units = c("in"), dpi = 300)
# If we want to see each state separately
jpeg(filename = "E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/trajectory/psedotime_heatmap.jpg",
     width=7,height = 9 ,units="in",res=1200)
plot_cell_trajectory(HSMM2, color_by = "Pseudotime",cell_size = 1) +
  facet_wrap(~State, nrow = 1)

plot_cell_trajectory(HSMM2, color_by = "State",cell_size = 1)+
  
  facet_wrap(~Pseudotime, nrow = 1)
plot_cell_trajectory(HSMM2, color_by = "Source",cell_size = 1,nrow=3)+
  
  facet_wrap(~Pseudotime, nrow = 1)

# DifferentialGeneTest for pseudotime

to_be_tested <- row.names(subset(fData(HSMM2),
                                 gene_short_name %in% obj.all.four_marker$X))
cds_subset <- HSMM2[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval <0.01)
sig_genes=read.csv("E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/trajectory/final/sigGenecluster_135.csv")
dim(sig_genes)
dim(sig_genes[,c("gene_short_name", "pval", "qval")])

#------------Pseudotemporal Expression Pattern
sig_gene_names <- row.names(subset(diff_test_res, qval <0.01))
#write.csv(sig_genes, "E:/thesis/MSc/final_sig_gene_trajectory.csv")
View(sig_gene_names)



k=plot_pseudotime_heatmap(HSMM2[sig_gene_names,],
                        num_clusters = 5,
                        cores = 1,
                        show_rownames = F)

five.cluster=c("TCL1A","IGHD")
three.cluster=c("CREM","DNAAF1","RGS1","ANKRD28","ANKRD37","IGHG1","IGKC","PPP1R15A","RRBP1","IGHGP","IGHG3","HSPD1","HSPH1","HSP90AA1","HSPA6","HSPB1","HSPA1A","HSPA1B","DNAJB1","IGHG4","IGKV1-5","MT1X")
two.cluster=c("IFRD1","C1orf56","CD69","CDC42SE1","LTB","HLA-DRB5","COTL1","CORO1A","IFITM2","RPS26","UQCRB","NBEAL1","ZNF90","AC090498.1","TRAC")

#jpeg(filename = "E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/trajectory/final/final_psedotime_heatmap_branch.jpg",
     width=4,height = 4 ,units="in",res=1200)
t=plot_genes_branched_heatmap(HSMM2[sig_gene_names,],
                            
                            num_clusters = 5,
                            cores = 1,
                            use_gene_short_name = F,
                            show_rownames = F)
#dev.cur()

t=plot_genes_branched_heatmap(HSMM2[sig_gene_names,],
                              
                              num_clusters = 5,
                              cores = 1,
                             return_heatmap = T)

 library(dplyr)
cluster_1_3_5=filter((t$annotation_row), Cluster%in%c(1,3,5))
[1] "FCN3"        "IFI6"        "JUN"         "GADD45A"     "MCL1"        "HSPA6"       "RGS1"        "ATF3"        "GUK1"        "ERLEC1"     
[11] "IGKC"        "IGKV1-12"    "IGKV4-1"     "IGKV3-11"    "IGKV3-20"    "MIR4435-2HG" "COBLL1"      "HSPD1"       "HSPE1"       "ANKRD28"    
[21] "SELK"        "MTRNR2L12"   "SSR3"        "JCHAIN"      "ANKRD37"     "MZB1"        "SQSTM1"      "HSPA1A"      "HSPA1B"      "CDKN1A"     
[31] "PRDM1"       "ZFAND2A"     "SEC61G"      "HSPB1"       "DNAJB9"      "PRDX4"       "PIM2"        "MTRNR2L10"   "NHSL2"       "SSR4"       
[41] "PTP4A3"      "POLR2L"      "MTRNR2L8"    "FKBP2"       "RAB30"       "ACTA2"       "DUSP5"       "FKBP11"      "TMED2"       "HSPH1"      
[51] "GAS6"        "FOS"         "HSP90AA1"    "IGHA2"       "IGHG4"       "IGHG2"       "IGHGP"       "IGHA1"       "IGHG1"       "IGHG3"      
[61] "IGHM"        "VIMP"        "MT2A"        "MT1E"        "MT1G"        "MT1X"        "DNAAF1"      "RRBP1"       "DNAJB1"      "PPP1R15A"   
[71] "SDF2L1"      "IGLV3-1"     "IGLC2"       "IGLC6"       "IGLC7"       "IGLL1"       "DERL3"       "XBP1"        "SELM"        "IGKV1-5"    
[81] "IGHV1-46" 
marker_orig.ident_cluster_1_2_4=filter(marker_orig.ident,gene%in%cluster_1_2_4)
dim(marker_orig.ident_cluster_1_2_4)
dev.off()

t$
t1 <- as.data.frame(cutree(k$tree_row, k=5))
                   colnames(t) <- "Cluster"
                   t$Gene <- rownames(t)
                 
df_cluster_135=filter(sig_genes, rownames(sig_genes)%in%rownames(cluster_1_3_5))
vargene_tested <- row.names(subset(fData(HSMM2),
gene_short_name %in% marker_orig.ident$gene))
plot_genes_in_pseudotime(HSMM2[c("IGKC","JCHAIN","IGHG1","IGHG2", "IGHG3","IGHG4"),],ncol=2,color_by = "Source")
selected=c("IGKV4-1","IGHG3","IGHG4","IGLV3-1","IGKC")

plot_genes_branched_pseudotime(HSMM2[c("IGHG3","IGHG4","IGHGP","IGHD","TRAC"),],
                               branch_point = 2,
                               color_by = "Source",
                               ncol = 1)
ggsave("./trajectory/final/selected_final_trajectory2col.jpg", width = 8, height = 4.8, units = c("in"), dpi = 300)

row.names(subset(fData(HSMM2),
                 gene_short_name %in% marker_orig.ident$gene))
marker_orig.ident = subset(marker_orig.ident,gene%in%sig_gene_names) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)


#===============pathway analysis

library(dplyr)
library(stringr)
library(KEGG.db)
library(fgsea)
library(GO.db)
library(ggfittext)
library(Seurat)
library(clusterProfiler)
library(purrr)


library(org.Hs.eg.db)
# get gene symbols; human = 'Hs'
SYMBOL2EG <-
  eval(parse(text = sprintf(
    'org.%s.egSYMBOL2EG', 'Hs'
  )))
#two.cluster_data=rownames(cluster_1_3_5)
two.cluster_data=filter(obj.all.four_marker,X %in%rownames(cluster_1_3_5))
marker_orig.ident=read.csv("orig.ident_marker.csv")
marker_orig.ident = marker_orig.ident %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)
View(obj.all.four@meta.data)
two.cluster_data=marker_orig.ident
logFC_score <- two.cluster_data$avg_log2FC
logFC_score
de_genes<- (two.cluster_data)$gene
length(x = de_genes)
names(logFC_score) <- de_genes
logFC_score

genes <- intersect(de_genes, mappedkeys(SYMBOL2EG)) #  Get the gene symbol that are mapped to an entrez gene identifiers
#length(x = genes)


logFC_score<- logFC_score[genes] # access logFC score 
#length(x = logFC_score) 
gene_entrez <-genes %>% SYMBOL2EG[.] %>% as.list %>% map( ~ .[1]) %>% simplify
#length(x = gene_entrez)

names(logFC_score) <- gene_entrez
logFC_score # input for rank in fgsea

#===========================GO pathway retrieving =========

library(gage)

go.hs=go.gsets(species='human')
go.bp=go.hs$go.sets[go.hs$go.subs$BP] # Biological Process
#length(x = go.bp) # 15,927
go.bp

# .... convert entrez gene id to gene Symbol
gte<- names(go.bp)
txp<- substr(gte, 12, nchar(gte))
substr(txp, 1, 1) <- toupper(substr(txp, 1, 1))
names(go.bp)<- txp
go.bp


fgseaRes.go <- fgseaMultilevel(go.bp, logFC_score,
                               scoreType = "pos",
                               eps      = 0.0,
                               minSize=1,
                               maxSize=Inf,
                               nPermSimple = 10000) # run fgsea

fgseaRes.go <- fgsea(go.bp, logFC_score,
                               scoreType = "pos",
                               minSize=1,
                               eps=1e-05,
                                maxSize=Inf,
                               nPerm = 10000)
#view(fgseaRes.go)
fgseaRes.go[, 7:8]
1e-05
# ===== plotting Upregulated pathways go

gu_ggdat_t <- fgseaRes.go %>% as.data.frame %>% as_data_frame %>% arrange((padj))
gu_ftr_res<- dplyr::filter(gu_ggdat_t,padj<=0.005)

library(ggplot2)
library(tidyverse)
# select top 10 paths
gu_ggdat<- gu_ftr_res %>% head(30) %>% mutate(pathway = fct_inorder(pathway))
padj_log10_kd<- round(-log(gu_ggdat$padj, base = 10),digits = 2)
# Pathway plot
ggplot(gu_ggdat) +
  geom_point(aes(
    x = (padj),
    y = pathway,
    #size = size,
    colour = size
  ), size = 8) + #geom_text(aes(x = pathway, y = NES, label = size, colour = size),
  # vjust = "inward", hjust = "inward",
  # show.legend = FALSE, size = 5, angle = 45)+
  labs(title = paste("GO Pathways", ": Upregulated"),
       x = '-log10(padj)',
       y = 'Pathways') +
  #scale_size_continuous(name = 'Size of\nthe pathway',range = c(2,8)) +
  theme_grey(base_size =14 ) +
  theme(axis.text.x = element_text(angle = -23, hjust = 0, size = 11, color  = "Black")
        ,plot.margin = margin(10,10,5,5)
  ) + scale_color_gradient("Size of\nthe pathway" , low="blue", high="red")


dotplot(gu_ggdat, showCategory=10) + facet_grid(.~.sign)

library(clusterProfiler)
library(KEGG.db)
library(plyr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggrepel)
library(fgsea)

library(reactome.db)
library(org.Hs.eg.db)
library(data.table)
require(DOSE)
library(DOSE)
library(clusterProfiler)
library(enrichplot)
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(org.Hs.eg.db)

gse <- gseGO(geneList=sort(logFC_score, decreasing = T), 
             ont ="BP", 
             
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db
          )
pathway_id=c("GO:0006911","GO:0006910","GO:0002455","GO:0006959","GO:0050853","GO:0019724","GO:0050864","GO:0042113",
  "GO:0002250","GO:0016064","GO:0002449","GO:0002682","GO:0002764","GO:0002253","GO:0002429","GO:0002768","GO:0042742","GO:0006956","GO:0006952","GO:0006955")
gse=gse$
gse1=(dplyr::filter(gse,ID%in% pathway_id))
View(gse@result)
#saveRDS(gse,"E:/thesis/MSc/msc_figure/bcell_from_all_datasets/trajectory/Final_gseGo_obj_for_trajectory_cluster_135.rds")
gse=readRDS("E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/trajectory/final/Final_gseGo_obj_for_trajectory_cluster_135.rds")
#write.csv(gse@result,"E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/trajectory/Final_pathway_gseGo_for_trajectory_cluster_135.csv")
gse=gse%>%arrange(desc(enrichmentScore))
p=DOSE::dotplot(gse, showCategory=25) 
df=read.csv("E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/trajectory/final/Final_pathway_gseGo_for_trajectory_cluster_135.csv")
df$GeneRatio=df$enrichmentScore
df$p.adjust=-log10(df$p.adjust)
#options(scipen = 999)
typeof(df$p.adjust)

library(ggplot2)
library(stringr)

ggplot(df,aes(
  x=reorder(Description,GeneRatio), 
  y=GeneRatio,
  color = p.adjust, size = GeneRatio
)) +
  geom_point()+
  scale_color_gradient(low = "red", high = "blue")+coord_flip()+theme_bw()+

  scale_x_discrete(labels = function(x) str_wrap(x, width =70))+
  theme(axis.text=element_text(color="black",size=13),
        panel.border = element_rect(color="black",fill=NA,linewidth =0.5),
        )+
  labs(title = "",
      x="",
       y = 'GeneRatio')

  

selected_pathway=gse[c("GO:0050853","GO:0002253","GO:0019724","GO:0002455","GO:0006959","GO:0050864","GO:0042113","GO:0006955")]
#write.csv(selected_pathway,"./pathway/selected_pathway_gsego_1st50.csv")
gene_pathway=(strsplit(selected_pathway$core_enrichment, "/"))
gene_pathway_list=c("3503", "3502", "3494", "3514", "3493", "3507", "3543", "3542", "3500", 
                    "3503",  "3502",  "3494",  "3320" , "3514",  "3493",  "3507",  "28299", "3543",  "3542"  ,"3500" ,
                    "3503",  "3502",  "3494",  "3514" , "3493",  "3507" , "28299", "3543" , "3542",  "3329" , "3500" , 
                    "3503",  "3502",  "3494",  "3514" , "3493",  "3507",  "28299", "3543",  "3542",  "3500" ,
                    "3503", "3502" ,"3494", "3514", "3493", "3507", "3543", "3542", "3500","3512",
                    "3503", "3502" ,"3494" ,"3514", "3493" ,"3507", "3543", "3542", "3329" ,"3500",
                    "3503",  "3502",  "3494" , "3303",  "3320" , "3304" , "3514" , "3310" , "3493" , "3507" , "28299" ,"3543"  ,"3542" , "3329" )
#ggsave("./trajectory/final/final_dotplot_goPathway_cluster123_with_25.jpg", width = 9, height = 7, units = c("in"), dpi = 300)

gene_pathway_unique=unique(gene_pathway_list)
mapped_Df=as.data.frame(mapIds(org.Hs.eg.db, c(gene_pathway_unique),"SYMBOL","ENTREZID"))
View(mapped_Df)
colnames(mapped_Df)="gene"
mapped_Df <- tibble::rownames_to_column(mapped_Df, "Entrez Id")
mapped_Df= read.csv(("./pathway/all_selected_pathway_gene_1st50.csv"))
#write.csv(mapped_Df ,"./pathway/all_selected_pathway_gene_1st50.csv")
library(dplyr)
merged=merge(mapped_Df,marker_orig.ident[,c("gene","avg_log2FC","p_val_adj","cluster")],by="gene")
merged=read.csv("./pathway/all_selected_pathway_merged_Df.csv")
#write.csv(merged ,"all_selected_pathway_merged_Df.csv")
selected=c("IGKV4-1","IGHG3","IGHG4","IGLV3-1","IGKC")
sig.genes
diff_sig_From_trajectory=filter(sig_genes, rownames(sig_genes)%in%merged$gene)%>%arrange(pval)
#write.csv(diff_sig_From_trajectory ,"./pathway/diff_sig_From_trajectory_1st50.csv")


#survival analysis
library(survival)
?lung
as_tibble(lung)
lung <- as_tibble(lung)
lung
Surv
s <- Surv(lung$time, lung$status)
survfit(s~1)
survfit(Surv(time, status)~1, data=lung)
sfit <- survfit(Surv(time, status)~1, data=lung)
sfit <- survfit(Surv(time, status)~sex, data=lung)
sfit <- survfit(Surv(time, status)~sex, data=lung)
plot(sfit)
