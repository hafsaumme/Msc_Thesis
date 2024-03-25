
library(stringr)

library(sctransform)
library(glmGamPoi)


#==========preparing GSE125449.set1_2_HCC data object=============

colnames(GSE125449.set1_2_HCC) <- sub('(.*?)_(.*)', '\\2_\\1', colnames(GSE125449.set1_2_HCC))

obj.GSE125449.set1_2_HCC = CreateSeuratObject(counts = GSE125449.set1_2_HCC, project ="GSE125449.set1_2_HCC", 
                                              min.cells = 3,
                                              min.features= 300)

obj.GSE125449.set1_2_HCC@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_HCC@meta.data), "H18"))] <- "H18"
obj.GSE125449.set1_2_HCC@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_HCC@meta.data), "H21"))] <- "H21"
obj.GSE125449.set1_2_HCC@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_HCC@meta.data), "H23"))] <- "H23"

obj.GSE125449.set1_2_HCC@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_HCC@meta.data), "H28"))] <- "H28"

obj.GSE125449.set1_2_HCC@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_HCC@meta.data), "H30"))] <- "H30"

obj.GSE125449.set1_2_HCC@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_HCC@meta.data), "H37"))] <- "H37"
obj.GSE125449.set1_2_HCC@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_HCC@meta.data), "H38"))] <- "H38"

obj.GSE125449.set1_2_HCC@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_HCC@meta.data), "H34"))] <- "H34"
obj.GSE125449.set1_2_HCC@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_HCC@meta.data), "H65"))] <- "H65"

obj.GSE125449.set1_2_HCC = PercentageFeatureSet(obj.GSE125449.set1_2_HCC, pattern = "^MT-", col.name = "percent.mt")

obj.GSE125449.set1_2_HCC@meta.data$orig.ident <- 'HCC'
obj.GSE125449.set1_2_HCC <- SetIdent(obj.GSE125449.set1_2_HCC, value = 'HCC')


# Visualize QC metrics as a violin plot
VlnPlot(obj.GSE125449.set1_2_HCC, features = c("nFeature_RNA",
                                  "nCount_RNA", "percent.mt"), ncol = 3)&
                                  theme(plot.title = element_text(color="black",size=13))
                                  
 
obj.GSE125449.set1_2_HCC <- subset(obj.GSE125449.set1_2_HCC, 
                                   subset = nCount_RNA>700 &nFeature_RNA > 200 & percent.mt < 20)

obj.GSE125449.set1_2_HCC 

obj.GSE125449.set1_2_HCC.plot1 <- FeatureScatter(obj.GSE125449.set1_2_HCC , feature1 = "nCount_RNA", 
                                                 feature2 = "percent.mt")
obj.GSE125449.set1_2_HCC.plot2 <- FeatureScatter(obj.GSE125449.set1_2_HCC , feature1 = "nCount_RNA", 
                                                 feature2 = "nFeature_RNA")
print(obj.GSE125449.set1_2_HCC.plot1 + obj.GSE125449.set1_2_HCC.plot2)

obj.GSE125449.set1_2_HCC <- SCTransform(obj.GSE125449.set1_2_HCC, method = "glmGamPoi", 
                                   vars.to.regress = "percent.mt", variable.features.n = 3000,
                                   verbose = FALSE)

dim(obj.GSE125449.set1_2_HCC) 

# Identify the 10 most highly variable genes

top_10.GSE125449.set1_2_HCC <- head(VariableFeatures(obj.GSE125449.set1_2_HCC), 10)

# plot variable features with and without labels

GSE125449.set1_2_HCC.plot.1 <- VariableFeaturePlot(obj.GSE125449.set1_2_HCC)
GSE125449.set1_2_HCC.plot.2 <- LabelPoints(plot = GSE125449.set1_2_HCC.plot.1, points = top_10.GSE125449.set1_2_HCC, repel = TRUE)
GSE125449.set1_2_HCC.plot.1 + GSE125449.set1_2_HCC.plot.2


# cell type annotation

set1.samples<-read.table("materials/GSE125449/set1/GSE125449_Set1_samples.txt",
                         header=TRUE,
                         sep ='\t')
dim(set1.samples)
colnames(set1.samples)[2]= "Cell.Barcode"

set2.samples<-read.table("materials/GSE125449/set2/GSE125449_Set2_samples.txt",
                         header=TRUE,
                         sep ='\t') 

set1_set2_samples= rbind(set1.samples, set2.samples)%>% select(-Sample)
dim(set1_set2_samples)
set1_set2_samples$Cell.Barcode= gsub("-", ".", set1_set2_samples$Cell.Barcode)

hcc_metData= obj.GSE125449.set1_2_HCC@meta.data
rownames(hcc_metData)= gsub("\\_.*","",rownames(hcc_metData))
hcc_metData <- rownames_to_column(hcc_metData, "Cell.Barcode")%>% select(liver, Cell.Barcode)
dim(hcc_metData)

hcc_metData= left_join(hcc_metData,set1_set2_samples,By="Cell.Barcode")


obj.GSE125449.set1_2_HCC@meta.data$newCellName= gsub("\\_.*","",rownames(obj.GSE125449.set1_2_HCC@meta.data))
hcc_metData <- hcc_metData[match((obj.GSE125449.set1_2_HCC@meta.data$newCellName), hcc_metData$Cell.Barcode), ] 
dim(hcc_metData)


obj.GSE125449.set1_2_HCC@meta.data$Type<-hcc_metData$Type





#==========preparing GSE125449.set1_2_iCCA data object=============

colnames(GSE125449.set1_2_iCCA) <- sub('(.*?)_(.*)', '\\2_\\1', colnames(GSE125449.set1_2_iCCA))

obj.GSE125449.set1_2_iCCA = CreateSeuratObject(counts = GSE125449.set1_2_iCCA, project ="GSE125449.set1_2_iCCA", 
                                               min.cells = 3,
                                               min.features= 300)

obj.GSE125449.set1_2_iCCA@meta.data$liver<-NA
obj.GSE125449.set1_2_iCCA@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_iCCA@meta.data), "C25"))] <- "C25"
obj.GSE125449.set1_2_iCCA@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_iCCA@meta.data), "C26"))] <- "C26"
obj.GSE125449.set1_2_iCCA@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_iCCA@meta.data), "C29"))] <- "C29"
obj.GSE125449.set1_2_iCCA@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_iCCA@meta.data), "C35"))] <- "C35"
obj.GSE125449.set1_2_iCCA@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_iCCA@meta.data), "C39"))] <- "C39"
obj.GSE125449.set1_2_iCCA@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_iCCA@meta.data), "C42"))] <- "C42"
obj.GSE125449.set1_2_iCCA@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_iCCA@meta.data), "C46"))] <- "C46"
obj.GSE125449.set1_2_iCCA@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_iCCA@meta.data), "C56"))] <- "C56"
obj.GSE125449.set1_2_iCCA@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_iCCA@meta.data), "C60"))] <- "C60"
obj.GSE125449.set1_2_iCCA@meta.data$liver[which(str_detect(rownames(obj.GSE125449.set1_2_iCCA@meta.data), "C66"))] <- "C66"

obj.GSE125449.set1_2_iCCA = PercentageFeatureSet(obj.GSE125449.set1_2_iCCA, pattern = "^MT-", col.name = "percent.mt")

obj.GSE125449.set1_2_iCCA$orig.ident="iCCA"

obj.GSE125449.set1_2_iCCA <- SetIdent(obj.GSE125449.set1_2_iCCA, value = 'iCCA')

# Visualize QC metrics as a violin plot
VlnPlot(obj.GSE125449.set1_2_iCCA, features = c("nFeature_RNA",
                                                          "nCount_RNA", "percent.mt"), ncol = 3)&
                                        theme(plot.title = element_text(color="black",size=13))


obj.GSE125449.set1_2_iCCA <- subset(obj.GSE125449.set1_2_iCCA, 
                                    subset = nCount_RNA>700 &nFeature_RNA > 200 & percent.mt < 20)

dim(obj.GSE125449.set1_2_iCCA)

obj.GSE125449.set1_2_iCCA.plot1 <- FeatureScatter(obj.GSE125449.set1_2_iCCA , feature1 = "nCount_RNA", 
                                                  feature2 = "percent.mt")
obj.GSE125449.set1_2_iCCA.plot2 <- FeatureScatter(obj.GSE125449.set1_2_iCCA , feature1 = "nCount_RNA", 
                                                  feature2 = "nFeature_RNA")
print(obj.GSE125449.set1_2_iCCA.plot1 + obj.GSE125449.set1_2_iCCA.plot2)

obj.GSE125449.set1_2_iCCA <- SCTransform(obj.GSE125449.set1_2_iCCA, method = "glmGamPoi", 
                                         vars.to.regress = "percent.mt", variable.features.n = 3000,
                                         verbose = FALSE)

dim(obj.GSE125449.set1_2_iCCA) 

colnames(obj.GSE125449.set1_2_HCC)

# cell type annotation===============

set1.samples<-read.table("materials/GSE125449/set1/GSE125449_Set1_samples.txt",
                         header=TRUE,
                         sep ='\t')
dim(set1.samples)
colnames(set1.samples)[2]= "Cell.Barcode"

set2.samples<-read.table("materials/GSE125449/set2/GSE125449_Set2_samples.txt",
                         header=TRUE,
                         sep ='\t') 

set1_set2_samples= rbind(set1.samples, set2.samples)%>% select(-Sample)
dim(set1_set2_samples)
set1_set2_samples$Cell.Barcode= gsub("-", ".", set1_set2_samples$Cell.Barcode)

iCCA_metData= obj.GSE125449.set1_2_iCCA@meta.data
rownames(iCCA_metData)= gsub("\\_.*","",rownames(iCCA_metData))
iCCA_metData <- rownames_to_column(iCCA_metData, "Cell.Barcode")%>% select(liver, Cell.Barcode)

dim(iCCA_metData)

iCCA_metData= left_join(iCCA_metData,set1_set2_samples,By="Cell.Barcode")

obj.GSE125449.set1_2_iCCA@meta.data$newCellName= gsub("\\_.*","",rownames(obj.GSE125449.set1_2_iCCA@meta.data))
iCCA_metData <- iCCA_metData[match((obj.GSE125449.set1_2_iCCA@meta.data$newCellName), iCCA_metData$Cell.Barcode), ] 
dim(iCCA_metData)


obj.GSE125449.set1_2_iCCA@meta.data$Type<-iCCA_metData$Type



#....................INTEGRATION on HCC and iCCA datasets...............

list= c(obj.GSE125449.set1_2_HCC, obj.GSE125449.set1_2_iCCA)


# select features that are repeatedly variable across datasets for integration

features <- SelectIntegrationFeatures(object.list = list, nfeatures = 4000)

list <- PrepSCTIntegration(object.list = list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features)
obj.hcc_icca <- IntegrateData(anchorset = anchors,normalization.method = "SCT")




obj.hcc_icca=RunPCA(obj.hcc_icca, features = VariableFeatures(object = obj.hcc_icca))

ElbowPlot(obj.hcc_icca, ndims = 50)

pct <- obj.hcc_icca[["pca"]]@stdev / sum(obj.hcc_icca[["pca"]]@stdev) * 100

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



obj.hcc_icca=RunUMAP(obj.hcc_icca, reduction = "pca", dims = 1:pcs)

obj.hcc_icca@meta.data$orig.ident[which(str_detect(obj.hcc_icca@meta.data$orig.ident,"HCC"))] <- 'HCC'
obj.hcc_icca@meta.data$orig.ident[which(str_detect(obj.hcc_icca@meta.data$orig.ident,"iCCA"))] <- 'iCCA'

DimPlot(obj.hcc_icca, reduction = "umap", group.by = "liver")+theme(axis.title = element_text(size=13))+ggtitle("")

DimPlot(obj.hcc_icca, reduction = "umap", group.by = "Type")+theme(axis.title = element_text(size=13))+ggtitle("")

DefaultAssay(obj.hcc_icca) <- "RNA"


# Retrieving B cells from HCC and iCCA datasets

obj.hcc_icc.bcell<- subset(obj.hcc_icca, subset = Type == 'B cell')

obj.bcell.hcc= subset(obj.hcc_icc.bcell, subset = orig.ident == 'HCC')

obj.bcell.icca= subset(obj.hcc_icc.bcell, subset = orig.ident == 'iCCA')


