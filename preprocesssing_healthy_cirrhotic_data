library(Seurat)
library(dplyr)
library(patchwork)
library(Seurat)
library(stringr)
library(HGNChelper)
library(ScType)

#===================== Reading Healthy data ==========================#

healthy_1.data_dir <- 'GSM4041150_healthy1_cd45+'
list.files(healthy_1.data_dir)

healthy_1.data <- Read10X(data.dir = healthy_1.data_dir)
dim(healthy_1.data)
colnames(healthy_1.data)

healthy_2.data_dir <- 'GSM4041151_healthy1_cd45-A'
list.files(healthy_2.data_dir)

healthy_2.data <- Read10X(data.dir = healthy_2.data_dir)
dim(healthy_2.data)

healthy_3.data_dir <- 'GSM4041152_healthy1_cd45-B'
list.files(healthy_3.data_dir)

healthy_3.data <- Read10X(data.dir = healthy_3.data_dir)
dim(healthy_3.data)

healthy_4.data_dir <- 'GSM4041153_healthy2_cd45+'
list.files(healthy_4.data_dir)

healthy_4.data <- Read10X(data.dir = healthy_4.data_dir)
dim(healthy_4.data) #33694  6571

healthy_5.data_dir <- 'GSM4041154_healthy2_cd45-'
list.files(healthy_5.data_dir)

healthy_5.data <- Read10X(data.dir = healthy_5.data_dir)
dim(healthy_5.data) #33694  4688

healthy_6.data_dir <- 'GSM4041155_healthy3_cd45+'
list.files(healthy_6.data_dir)

healthy_6.data <- Read10X(data.dir = healthy_6.data_dir)
dim(healthy_6.data) #33694  3295

healthy_7.data_dir <- 'GSM4041156_healthy3_cd45-A'
list.files(healthy_7.data_dir)

healthy_7.data <- Read10X(data.dir = healthy_7.data_dir)
dim(healthy_7.data) #33694  3230

healthy_8.data_dir <- 'GSM4041157_healthy3_cd45-B'
list.files(healthy_8.data_dir)

healthy_8.data <- Read10X(data.dir = healthy_8.data_dir)
dim(healthy_8.data) #33694  1156

healthy_9.data_dir <- 'mGSM4041158_healthy4_cd45+'
list.files(healthy_9.data_dir)

healthy_9.data <- Read10X(data.dir = healthy_9.data_dir)
dim(healthy_9.data) #33694  4937

healthy_10.data_dir <- 'GSM4041159_healthy4_cd45-'
list.files(healthy_10.data_dir)

healthy_10.data <- Read10X(data.dir = healthy_10.data_dir)
dim(healthy_10.data) #33694  3306

healthy_11.data_dir <- 'GSM4041160_healthy5_cd45+'
list.files(healthy_11.data_dir)

healthy_11.data <- Read10X(data.dir = healthy_11.data_dir)
dim(healthy_11.data) #33694  5122


#============================= Reading Cirrhotic data=========================

cirrhotic_1.data_dir <- 'GSM4041161_cirrhotic1_cd45+'
list.files(cirrhotic_1.data_dir)

cirrhotic_1.data <- Read10X(data.dir = cirrhotic_1.data_dir)
dim(cirrhotic_1.data) #33694  2380

cirrhotic_2.data_dir <- 'GSM4041162_cirrhotic1_cd45-A'
list.files(cirrhotic_2.data_dir)

cirrhotic_2.data <- Read10X(data.dir = cirrhotic_2.data_dir)
dim(cirrhotic_2.data) #33694  1686

cirrhotic_3.data_dir <- 'GSM4041163_cirrhotic1_cd45-B'
list.files(cirrhotic_3.data_dir)

cirrhotic_3.data <- Read10X(data.dir = cirrhotic_3.data_dir)
dim(cirrhotic_3.data) #33694  2046

cirrhotic_4.data_dir <- 'GSM4041164_cirrhotic2_cd45'
list.files(cirrhotic_4.data_dir)

cirrhotic_4.data <- Read10X(data.dir = cirrhotic_4.data_dir)
dim(cirrhotic_4.data) #33694  3412

cirrhotic_5.data_dir <- 'GSM4041165_cirrhotic2_cd45-'
list.files(cirrhotic_5.data_dir)

cirrhotic_5.data <- Read10X(data.dir = cirrhotic_5.data_dir)
dim(cirrhotic_5.data) #33694  3673

cirrhotic_6.data_dir <- 'GSM4041166_cirrhotic3_cd45+'
list.files(cirrhotic_6.data_dir)

cirrhotic_6.data <- Read10X(data.dir = cirrhotic_6.data_dir)
dim(cirrhotic_6.data) #33694  1985

cirrhotic_7.data_dir <- 'GSM4041167_cirrhotic3_cd45-'
list.files(cirrhotic_7.data_dir)

cirrhotic_7.data <- Read10X(data.dir = cirrhotic_7.data_dir)
dim(cirrhotic_7.data) #33694  3602

cirrhotic_8.data_dir <- 'GSM4041168_cirrhotic4_cd45+'
list.files(cirrhotic_8.data_dir)

cirrhotic_8.data <- Read10X(data.dir = cirrhotic_8.data_dir)
dim(cirrhotic_8.data) #33694  4819

cirrhotic_9.data_dir <- 'GSM4041169_cirrhotic5_cd45+'
list.files(cirrhotic_9.data_dir)

cirrhotic_9.data <- Read10X(data.dir = cirrhotic_9.data_dir)
dim(cirrhotic_9.data) #33694  2922




#===========================Preprocessing ===========================================================


preprocessing <- function(dataCount,proj_test){
  
  data <- CreateSeuratObject(counts = dataCount, 
                             project =proj_test)
  data$Liver <- proj_test
  data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
  vln<-VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  print(vln)
  data.plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  
  data.plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  print(data.plot1+ data.plot2)
  
  data <- subset(data,  percent.mt < 30)
  data <- SCTransform(data, method = "glmGamPoi", return.only.var.genes=F,
                      vars.to.regress = "percent.mt",
                      verbose = FALSE)
  
  return(data)
}

healthy_1.obj <-preprocessing(healthy_1.data, "healthy_1")



healthy_2.obj <-preprocessing(healthy_2.data, "healthy_1")


healthy_3.obj <-preprocessing(healthy_3.data, "healthy_1")



healthy_4.obj <-preprocessing(healthy_4.data, "healthy_2")


healthy_5.obj <-preprocessing(healthy_5.data, "healthy_2")


healthy_6.obj <-preprocessing(healthy_6.data, "healthy_3")


healthy_7.obj <-preprocessing(healthy_7.data, "healthy_3")



healthy_8.obj <-preprocessing(healthy_8.data, "healthy_3")


healthy_9.obj <-preprocessing(healthy_9.data, "healthy_4")



healthy_10.obj <-preprocessing(healthy_10.data, "healthy_4")



healthy_11.obj <-preprocessing(healthy_11.data, "healthy_5")



#....................INTEGRATION oh healthy datasets...............

healthy.list= c(healthy_1.obj,healthy_2.obj,
                healthy_3.obj,healthy_4.obj,
                healthy_5.obj,healthy_6.obj,
                healthy_7.obj,healthy_8.obj,
                healthy_9.obj,healthy_10.obj,
                healthy_11.obj)

# select features that are repeatedly variable across datasets for integration
healthy.features <- SelectIntegrationFeatures(object.list = healthy.list,nfeatures=3000)
healthy.list <- PrepSCTIntegration(object.list = healthy.list, anchor.features = healthy.features)

healthy.anchors <- FindIntegrationAnchors(object.list = healthy.list, anchor.features = healthy.features)
# this command creates an 'integrated' data assay
healthy_int <- IntegrateData(anchorset = healthy.anchors)
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(healthy_int) <- "SCT"



healthy_int <- SCTransform(healthy_int, method = "glmGamPoi", 
                           vars.to.regress = "percent.mt",variable.features.n = 3000,
                           verbose = FALSE)

top_10.healthy_int<- head(top_10.healthy_int<- head(VariableFeatures(healthy_int), 10))

# plot variable features with and without labels
healthy_int.plot.1 <- VariableFeaturePlot((healthy_int))
healthy_int.plot.2 <- LabelPoints(plot = healthy_int.plot.1, points = top_10.healthy_int, repel = TRUE)
healthy_int.plot.1 + healthy_int.plot.2

healthy_int.vln<-VlnPlot(healthy_int, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)&
  theme(plot.title = element_text(color="black",size=13))

healthy_int.plot1 <- FeatureScatter(healthy_int, feature1 = "nCount_RNA", feature2 = "percent.mt")

healthy_int.plot2 <- FeatureScatter(healthy_int, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

print(data.plot1+ data.plot2)

healthy_int  <- RunPCA(healthy_int,npcs = 50,
                       features = VariableFeatures(object = healthy_int))

DimPlot(healthy_int, reduction = "pca")

DimHeatmap(healthy_int, dims = 1:10, cells = 500, balanced = TRUE)

# Determining useful PCA number ===============================================

ElbowPlot(healthy_int, ndims = 35)

# Determine percent of variation associated with each PC
pct <- healthy_int[["pca"]]@stdev / sum(healthy_int[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
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



healthy_int <- RunUMAP(healthy_int, reduction = "pca", 
                       dims = 1:pcs)

#============patients cluster==========
healthy_int$orig.ident="healthy"
p=DimPlot(healthy_int, reduction = "umap", group.by = "Type")+
  theme(axis.title = element_text(size=13))+ggtitle("")+labs(color="")


#===========================Preprocessing Cirrotic datasets =============================


cirrhotic_1.obj <-preprocessing(cirrhotic_1.data, "cirrhotic_1")


cirrhotic_2.obj <-preprocessing(cirrhotic_2.data, "cirrhotic_1")

cirrhotic_3.obj <-preprocessing(cirrhotic_3.data, "cirrhotic_1")


cirrhotic_4.obj <-preprocessing(cirrhotic_4.data, "cirrhotic_2")


cirrhotic_5.obj <-preprocessing(cirrhotic_5.data, "cirrhotic_2")


cirrhotic_6.obj <-preprocessing(cirrhotic_6.data, "cirrhotic_3")


cirrhotic_7.obj <-preprocessing(cirrhotic_7.data, "cirrhotic_3")


cirrhotic_8.obj <-preprocessing(cirrhotic_8.data, "cirrhotic_4")


cirrhotic_9.obj <-preprocessing(cirrhotic_9.data, "cirrhotic_5")


#....................INTEGRATION oh cirrhotic datasets...............

cirrhotic.list= c(cirrhotic_1.obj,cirrhotic_2.obj,
                  cirrhotic_3.obj,cirrhotic_4.obj,
                  cirrhotic_5.obj,cirrhotic_6.obj,
                  cirrhotic_7.obj,cirrhotic_8.obj,
                  cirrhotic_9.obj)

# select features that are repeatedly variable across datasets for integration
cirrhotic.features <- SelectIntegrationFeatures(object.list = cirrhotic.list)

cirrhotic.anchors <- FindIntegrationAnchors(object.list = cirrhotic.list, anchor.features = cirrhotic.features)
# this command creates an 'integrated' data assay
cirrhotic_int <- IntegrateData(anchorset = cirrhotic.anchors)
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(cirrhotic_int) <- "integrated"

cirrhotic_int <- SCTransform(cirrhotic_int, method = "glmGamPoi", 
                             vars.to.regress = "percent.mt",variable.features.n = 3000,
                             verbose = FALSE)

cirrhotic_int.vln<-VlnPlot(cirrhotic_int, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)&
  theme(plot.title = element_text(color="black",size=13))
print(vln.iCCA)

cirrhotic_int.plot1 <- FeatureScatter(cirrhotic_int, feature1 = "nCount_RNA", feature2 = "percent.mt")

cirrhotic_int.plot2 <- FeatureScatter(cirrhotic_int, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

print(data.plot1+ data.plot2)

top_10.cirrhotic_int<- head(VariableFeatures(cirrhotic_int), 10)

# plot variable features with and without labels
cirrhotic_int.plot.1 <- VariableFeaturePlot((cirrhotic_int))
cirrhotic_int.plot.2 <- LabelPoints(plot = cirrhotic_int.plot.1, points = top_10.cirrhotic_int, repel = TRUE)
cirrhotic_int.plot.1 + cirrhotic_int.plot.2

cirrhotic_int  <- RunPCA(cirrhotic_int , npcs = 50,
                         features = VariableFeatures(object = cirrhotic_int))

DimPlot(cirrhotic_int, reduction = "pca")

DimHeatmap(cirrhotic_int, dims = 1:10, cells = 500, balanced = TRUE)

# Determining useful PCA number ===============================================

#DimHeatmap(data.combined, dims = 1:10, cells =500, balanced = TRUE)
ElbowPlot(cirrhotic_int, ndims = 35)

# Determine percent of variation associated with each PC
pct <- cirrhotic_int[["pca"]]@stdev / sum(cirrhotic_int[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
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



cirrhotic_int <- RunUMAP(cirrhotic_int, reduction = "pca", 
                         dims = 1:pcs)


#============patients cluster==========
p=DimPlot(cirrhotic_int, reduction = "umap", group.by = "Liver")+
  theme(axis.title = element_text(size=13))+ggtitle("")+labs(color="")


#....................INTEGRATION oh healthy and cirrhotic datasets...............

list= c(healthy_1.obj,healthy_2.obj,
        healthy_3.obj,healthy_4.obj,
        healthy_5.obj,healthy_6.obj,
        healthy_7.obj,healthy_8.obj,
        healthy_9.obj,healthy_10.obj,
        healthy_11.obj, cirrhotic_1.obj,cirrhotic_2.obj,
        cirrhotic_3.obj,cirrhotic_4.obj,
        cirrhotic_5.obj,cirrhotic_6.obj,
        cirrhotic_7.obj,cirrhotic_8.obj,
        cirrhotic_9.obj)


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 2000)

list <- PrepSCTIntegration(object.list = list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features)
combined <- IntegrateData(anchorset = anchors,features.to.integrate = features_integrate, normalization.method = "SCT")
levels(combined@active.ident)=c(levels(combined@active.ident),"healthy")
levels(combined@active.ident)=c(levels(combined@active.ident),"cirrhotic")
combined@active.ident[which(str_detect((combined@active.ident), "healthy"))] <- "healthy"
combined@active.ident[which(str_detect((combined@active.ident), "cirrhotic"))] <- "cirrhotic"
combined@meta.data$orig.ident[which(str_detect((combined@meta.data$orig.ident), "healthy"))] <- "healthy"
combined@meta.data$orig.ident[which(str_detect((combined@meta.data$orig.ident), "cirrhotic"))] <- "cirrhotic"

DefaultAssay(combined) <- "integrated"

combined=RunPCA(combined,  features = VariableFeatures(object = combined))

ElbowPlot(combined, ndims = 50)

pct <- combined[["pca"]]@stdev / sum(combined[["pca"]]@stdev) * 100

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



combined=RunUMAP(combined, reduction = "pca", dims = 1:15)
DimPlot(combined, reduction = "umap")

combined <- FindNeighbors(combined, reduction = "pca", dims = 1:15)

combined= FindClusters(combined,resolution = 0.5)

DimPlot(combined, reduction = "umap")

#finding cell type

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "ScTypeDB_full.xlsx"; tissue = "Liver"

gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = combined@assays$RNA@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = NULL) 
#View(es.max)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(combined@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(combined@meta.data[combined@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(combined@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

combined@meta.data$customclassif = NA
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  combined@meta.data$customclassif[combined@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

#View(combined@meta.data)

DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif') 


# Rename all identities
combined <- RenameIdents(object = combined, 
                         "0" = "Tcell",
                         "1" = "Tcell",
                         "2" = "Tcell",
                         "3" = "ILC",
                         "4" = "MP",
                         "5" = "Endothelia",
                         "6" = "ILC",
                         "7" = "ILC",
                         "8" = "Mesenchyme",
                         "9" = "Endothelia",
                         "10" = "Epithelia",
                         "11" = "MP",
                         "12" = "MP",
                         "13" = "ILC",
                         "14" = "Bcell",
                         "15" = "Plasma",
                         "16" = "Cycling",
                         "17" = "MP", 
                         "18" = "MP", 
                         "19" = "pDC", 
                         "20" = "Endothelia",
                         "21" = "Mast")


DimPlot(combined, reduction = "umap", group.by = "Liver")+theme(plot.title=element_blank(),legend.text = element_text(size=14))

DefaultAssay(combined) <- "RNA"

combined.markers$cell_type=NA
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
combined.markers = combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 1000, order_by = avg_log2FC)
combined.markers$cell_type[which(combined.markers$cluster==21)] <- 'Mast'
combined.markers$cell_type[which(combined.markers$cluster==19)] <- 'pDC'
combined.markers$cell_type[which(combined.markers$cluster==16)] <- 'Cycling'
combined.markers$cell_type[which(combined.markers$cluster==15)] <- 'Plasma'
combined.markers$cell_type[which(combined.markers$cluster==10)] <- 'Epithelia'
combined.markers$cell_type[which(combined.markers$cluster==8)] <- 'Mesenchyme'
combined.markers$cell_type[which(combined.markers$cluster==14)] <- 'Bcell'
combined.markers$cell_type[which(combined.markers$cluster==0|combined.markers$cluster==1|combined.markers$cluster==2)] <- 'Tcell'
combined.markers$cell_type[which(combined.markers$cluster==3|combined.markers$cluster==6|combined.markers$cluster==7|combined.markers$cluster==13)] <- 'ILC'
combined.markers$cell_type[which(combined.markers$cluster==4|combined.markers$cluster==11|combined.markers$cluster==12|combined.markers$cluster==17|combined.markers$cluster==18)] <- 'MP'
combined.markers$cell_type[which(combined.markers$cluster==5|combined.markers$cluster==9|combined.markers$cluster==20)] <- 'Endothelia'
combined.markers = combined.markers %>%
  group_by(cell_type) %>%
  slice_max(n = 1000, order_by = avg_log2FC)
combined@assays$RNA@data=NormalizeData(combined@assays$RNA@counts, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(combined)
combined@assays$RNA@scale.data= ScaleData(combined@assays$RNA@data, features = all.genes)




mp=c("FCN1","VCAN",
     "C1QC","IL1B","MNDA","CPVL","MS4A7","MSR1","CD163","CYBB",
     "CSF2RA","CLEC10A")

endothelia= c("CLEC14A","PLVAP","EGFL7","VWF","HSPG2","EMCN",
              "FCN3","OIT3","CLEC4G","CLEC4M")

mesenchyme=c("TAGLN","RGS5","ACTA2","PDGFRB",
             "TIMP1","LUM","COL3A1","PDGFRA","COL1A2","CCL19")

epithelia=c("ALB","F2",
            "KNG1","PRAP1","TTR","HSD17B6","ORM1","ITIH1","TF","CYP4A11",
            "APOC2","AGXT","KRT19","CXCL1","CXCL6","SOX9")
bcell=c("CD79A","MS4A1")
plasma=c("FCRL5","JSRP1")
cycling=c("MKI67","TOP2A","UBE2C","RRM2")
pDC=c("LILRA4","LRRC26","PTCRA","CLEC4C")
mast=c("TPSAB1","CPA3")
DefaultAssay(combined) <- "RNA"
subset=subset(x = combined, downsample = 1000)
all.genes <- rownames(subset)

DoHeatmap(subset,features = top10$gene) + NoLegend()
DoHeatmap(subset,features = marker,size=3,angle=40)+theme(axis.text=element_text(size=9,color="black"))

gene_p=c("TRAC","CD2","KLRF1","KLRC1","PRF1","FGFBP2","FCN1","VCAN",
         "C1QC","IL1B","MNDA","CPVL","MS4A7","MSR1","CD163","CYBB",
         "CSF2RA","CLEC10A","CLEC14A","PLVAP","EGFL7","VWF","HSPG2","EMCN",
         "FCN3","OIT3","CLEC4G","CLEC4M","TAGLN","RGS5","ACTA2","PDGFRB",
         "TIMP1","LUM","COL3A1","PDGFRA","COL1A2","CCL19","ALB","F2",
         "KNG1","PRAP1","TTR","HSD17B6","ORM1","ITIH1","TF","CYP4A11",
         "APOC2","AGXT","KRT19","CXCL1","CXCL6","SOX9","CD79A","MS4A1","FCRL5","JSRP1",
         "MKI67","TOP2A","UBE2C","RRM2","LILRA4","LRRC26","PTCRA","CLEC4C","TPSAB1","CPA3")
gene_paper=data.frame("gene"=gene_p)

gene_paper$cluster=NA
gene_paper$cluster[1:2]="Tcell"
gene_paper$cluster[3:6]="ILC"
gene_paper$cluster[7:18]="MP"
gene_paper$cluster[19:28]="Endothelia"
gene_paper$cluster[29:38]="Mesenchyme"
gene_paper$cluster[39:54]="Epithelia"
gene_paper$cluster[55:56]="Bcell"
gene_paper$cluster[57:58]="Plasma"
gene_paper$cluster[59:62]="Cycling"
gene_paper$cluster[63:66]="pDC"
gene_paper$cluster[67:68]="Mast"
library(ggplot2)
DoHeatmap(subset,features = top10$gene) + NoLegend()
DoHeatmap(subset,features = (gene_p))+theme(axis.text.y = element_text(size = 4))


healthy_cirrhotic_Bcell<- subset((combined), subset = customclassif == 'Bcell')
healthy_cirrhotic_Plasma<- subset((combined), subset = customclassif == 'Plasma')



