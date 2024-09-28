
library(Seurat)
library(ggplot2)
library(patchwork)
library(tibble)
library(dplyr)

GSE125449.set1 <- Read10X(data.dir = "materials/GSE125449/set1/") 
 
GSE125449.set1.samples<-read.table("materials/GSE125449/set1/GSE125449_Set1_samples.txt",
                                   header=TRUE,
                                   sep ='\t')


# ===========================filtering HCC and iCCA data from  GSE125449.set1.samples ==============================================

GSE125449.set1.samples.HCC<-filter(GSE125449.set1.samples, Sample %in% 
                                     c("S02_P01_LCP21","S10_P05_LCP23",
                                       "S07_P02_LCP28","S12_P07_LCP30",
                                       "S21_P13_LCP37", "S15_P09_LCP38","S16_P10_LCP18"))
dim(GSE125449.set1.samples.HCC) # 3086    3

#converting GSE125449.set1 to matrix

GSE125449.set1.frame= (as.data.frame(GSE125449.set1))

#finding subset for all HCC sample from set1

GSE125449.set1.HCC<- GSE125449.set1.frame[, GSE125449.set1.samples.HCC$Cell_Barcode]


#finding subset for each HCC cell_barcode from set2 

set1.HCC_sample= c("S02_P01_LCP21","S10_P05_LCP23",
                   "S07_P02_LCP28","S12_P07_LCP30",
                   "S21_P13_LCP37", "S15_P09_LCP38","S16_P10_LCP18")

set1.HCC_prefix_barcode=c("H21","H23","H28","H30","H37","H38","H18")

set1.HCC_variable=c("GSE125449.set1.H21","GSE125449.set1.H23",
                    "GSE125449.set1.H28","GSE125449.set1.H30",
                    "GSE125449.set1.H37","GSE125449.set1.H38",
                    "GSE125449.set1.H18")


j=1
for(samples in set1.HCC_sample) 
  
{
  
  filter_cell = assign(set1.HCC_variable[j], GSE125449.set1.frame[, filter(GSE125449.set1.samples, Sample %in% samples)$Cell_Barcode])
  print(samples)
  
  j=j+1
  
}
#adding prefix with cell_barcode HCC set1
colnames(GSE125449.set1.H21)= paste(set1.HCC_prefix_barcode[1], colnames(GSE125449.set1.H21), sep = "_")
colnames(GSE125449.set1.H23)= paste(set1.HCC_prefix_barcode[2], colnames(GSE125449.set1.H23), sep = "_")
colnames(GSE125449.set1.H28)= paste(set1.HCC_prefix_barcode[3], colnames(GSE125449.set1.H28), sep = "_")
colnames(GSE125449.set1.H30)= paste(set1.HCC_prefix_barcode[4], colnames(GSE125449.set1.H30), sep = "_")
colnames(GSE125449.set1.H37)= paste(set1.HCC_prefix_barcode[5], colnames(GSE125449.set1.H37), sep = "_")
colnames(GSE125449.set1.H38)= paste(set1.HCC_prefix_barcode[6], colnames(GSE125449.set1.H38), sep = "_")
colnames(GSE125449.set1.H18)= paste(set1.HCC_prefix_barcode[7], colnames(GSE125449.set1.H18), sep = "_")

dim(GSE125449.set1.H21) #20124   704
dim(GSE125449.set1.H23) #20124   151
dim(GSE125449.set1.H28) #20124   124
dim(GSE125449.set1.H30) #20124   805
dim(GSE125449.set1.H37) #20124   132
dim(GSE125449.set1.H38) #20124  1046
dim(GSE125449.set1.H18) #20124   124


GSE125449.set1.prefix_HCC=(data.frame(GSE125449.set1.H21, GSE125449.set1.H23,
                                      GSE125449.set1.H28, GSE125449.set1.H30,
                                      GSE125449.set1.H37, GSE125449.set1.H38,
                                      GSE125449.set1.H18))

dim((GSE125449.set1.prefix_HCC)) #20124  3086




#### end #####

GSE125449.set1.samples.iCCA<-filter(GSE125449.set1.samples, Sample %in% 
                                      c("S09_P04_LCP25","S08_P03_LCP26",
                                        "S11_P06_LCP29","S20_P12_LCP35",
                                        "S19_P11_LCP39"))

#finding subset all iCCA from set1
GSE125449.set1.iCCA<- GSE125449.set1.frame[, GSE125449.set1.samples.iCCA$Cell_Barcode]


#finding subset for each iCCA cell_barcode from set1 

set1.iCCA_sample= c("S09_P04_LCP25","S08_P03_LCP26",
                    "S11_P06_LCP29","S20_P12_LCP35",
                    "S19_P11_LCP39")

set1.iCCA_prefix_barcode=c("C25","C26","C29","C35","C39")

set1.iCCA_variable=c("GSE125449.set1.C25","GSE125449.set1.C26",
                     "GSE125449.set1.C29","GSE125449.set1.C35",
                     "GSE125449.set1.C39")

i=1

for(samples_c in set1.iCCA_sample) 
  
{
  
  filter_cell = assign(set1.iCCA_variable[i], GSE125449.set1.frame[, filter(GSE125449.set1.samples, Sample %in% samples_c)$Cell_Barcode])
  print(samples_c)
  i=i+1
  
}

#adding prefix with cell_barcode iCCA set1
colnames(GSE125449.set1.C25)= paste(set1.iCCA_prefix_barcode[1], colnames(GSE125449.set1.C25), sep = "_")
colnames(GSE125449.set1.C26)= paste(set1.iCCA_prefix_barcode[2], colnames(GSE125449.set1.C26), sep = "_")
colnames(GSE125449.set1.C29)= paste(set1.iCCA_prefix_barcode[3], colnames(GSE125449.set1.C29), sep = "_")
colnames(GSE125449.set1.C35)= paste(set1.iCCA_prefix_barcode[4], colnames(GSE125449.set1.C35), sep = "_")
colnames(GSE125449.set1.C39)= paste(set1.iCCA_prefix_barcode[5], colnames(GSE125449.set1.C39), sep = "_")


dim(GSE125449.set1.C25) #20124   207
dim(GSE125449.set1.C26) #20124   299
dim(GSE125449.set1.C29) #20124   939
dim(GSE125449.set1.C35) #20124   139
dim(GSE125449.set1.C39) #20124   445

GSE125449.set1.prefix_iCCA=(data.frame(GSE125449.set1.C25, GSE125449.set1.C26,
                                       GSE125449.set1.C29, GSE125449.set1.C35,
                                       GSE125449.set1.C39))
dim(GSE125449.set1.prefix_iCCA)

# ===========================filtering HCC and iCCA data from  GSE125449.set2.samples ==============================================

GSE125449.set2 <- Read10X(data.dir = "materials/GSE125449/set2/")

dim(GSE125449.set2) # 19572 X 4831

GSE125449.set2.samples<-read.table("materials/GSE125449/set2/GSE125449_Set2_samples.txt",
                                   header=TRUE,
                                   sep ='\t')

# filtering HCC data from  GSE125449.set2.samples----------------------------

GSE125449.set2.samples.HCC<-filter(GSE125449.set2.samples, Sample %in% c("S351_P10_LCP34","S364_P21_LCP65"))
dim(GSE125449.set2.samples.HCC) #827   3

#converting GSE125449.set2 to matrix

GSE125449.set2.frame= (as.data.frame(GSE125449.set2))

#finding subset for all HCC cell_barcode from set2

GSE125449.set2.HCC<- GSE125449.set2.frame[, GSE125449.set2.samples.HCC$Cell.Barcode]

dim(GSE125449.set2.HCC) #19572   827

#finding subset for each HCC cell_barcode from set2 

GSE125449.set2.H34<- GSE125449.set2.frame[, filter(GSE125449.set2.samples, Sample %in% c("S351_P10_LCP34"))$Cell.Barcode]

GSE125449.set2.H65<- GSE125449.set2.frame[, filter(GSE125449.set2.samples, Sample %in% c("S364_P21_LCP65"))$Cell.Barcode]

#adding suffix H34 with cell_barcode

colnames(GSE125449.set2.H34) <- paste("H34", colnames(GSE125449.set2.H34), sep = "_")

colnames(GSE125449.set2.H65) <- paste("H65", colnames(GSE125449.set2.H65), sep = "_")

GSE125449.set2.prefix_HCC= data.frame(GSE125449.set2.H34, GSE125449.set2.H65)

dim(GSE125449.set2.prefix_HCC)

# filtering iCCA data from  GSE125449.set2.samples----------------------------

GSE125449.set2.samples.iCCA<-filter(GSE125449.set2.samples, Sample %in% c("S355_P13_LCP42","S358_P16_LCP46",
                                                                          "S305_P06_LCP56","S300_P02_LCP60",
                                                                          "S365_P22_LCP66"))
dim(GSE125449.set2.samples.iCCA) #4004    3

#finding subset for all iCCA cell_barcode from set2

GSE125449.set2.iCCA<- GSE125449.set2.frame[, GSE125449.set2.samples.iCCA$Cell.Barcode]

length(GSE125449.set2.iCCA) #4004

#finding subset for each iCCA cell_barcode from set2 

GSE125449.set2.C42<- GSE125449.set2.frame[, filter(GSE125449.set2.samples, Sample %in% c("S355_P13_LCP42"))$Cell.Barcode]

#adding suffix C42 with cell_barcode

colnames(GSE125449.set2.C42) <- paste("C42", colnames(GSE125449.set2.C42), sep = "_")

GSE125449.set2.C46<- GSE125449.set2.frame[, filter(GSE125449.set2.samples, Sample %in% c("S358_P16_LCP46"))$Cell.Barcode]

length(GSE125449.set2.C46) #585

#adding suffix C46 with cell_barcode

colnames(GSE125449.set2.C46) <- paste("C46", colnames(GSE125449.set2.C46), sep = "_")

GSE125449.set2.C56<- GSE125449.set2.frame[, filter(GSE125449.set2.samples, Sample %in% c("S305_P06_LCP56"))$Cell.Barcode]

length( GSE125449.set2.C56) #137

#adding suffix C56 with cell_barcode
colnames( GSE125449.set2.C56) <- paste("C56", colnames(GSE125449.set2.C56), sep = "_")

GSE125449.set2.C60<- GSE125449.set2.frame[, filter(GSE125449.set2.samples, Sample %in% c("S300_P02_LCP60"))$Cell.Barcode]

length(GSE125449.set2.C60) #1418

#adding suffix C60 with cell_barcode

colnames(GSE125449.set2.C60) <- paste("C60", colnames(GSE125449.set2.C60), sep = "_")

GSE125449.set2.C66<- GSE125449.set2.frame[, filter(GSE125449.set2.samples, Sample %in% c("S365_P22_LCP66"))$Cell.Barcode]

length(GSE125449.set2.C66) #1356

#adding suffix C66 with cell_barcode

colnames(GSE125449.set2.C66) <- paste("C66", colnames(GSE125449.set2.C66), sep = "_")

GSE125449.set2.prefix_iCCA=(data.frame(GSE125449.set2.C42, GSE125449.set2.C46,
                                       GSE125449.set2.C56, GSE125449.set2.C60,
                                       GSE125449.set2.C66))
dim(GSE125449.set1.prefix_iCCA)

#===============Merging set1 and set2 for iCCA===========

common_gene= intersect(rownames(GSE125449.set1.prefix_iCCA), rownames(GSE125449.set2.prefix_iCCA))

a <- tibble::rownames_to_column(GSE125449.set2.prefix_iCCA, "row.names")

b<-  tibble::rownames_to_column(GSE125449.set1.prefix_iCCA, "row.names")

GSE125449.set1_2_iCCA= inner_join(a, b, by="row.names")

row.names(GSE125449.set1_2_iCCA) <- GSE125449.set1_2_iCCA$row.names

GSE125449.set1_2_iCCA[1] <- NULL


#================ Merging set1 and set2 for HCC===========


common_geneHCC= intersect(rownames(GSE125449.set1.prefix_HCC), rownames(GSE125449.set2.prefix_HCC))

a1 <- tibble::rownames_to_column(GSE125449.set2.prefix_HCC, "row.names")

b1<-  tibble::rownames_to_column(GSE125449.set1.prefix_HCC, "row.names")
 
GSE125449.set1_2_HCC= inner_join(a1, b1, by="row.names")

row.names(GSE125449.set1_2_HCC) <- GSE125449.set1_2_HCC$row.names

GSE125449.set1_2_HCC[1] <- NULL

dim(GSE125449.set1_2_HCC)


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
