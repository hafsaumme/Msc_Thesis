
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



