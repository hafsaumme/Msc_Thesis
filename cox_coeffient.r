library(readxl)
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(tibble)
library("stringr")
library(stringr)
library(glmGamPoi)

#write.csv(obj.all.four_marker,"materials/hafsa/uh_msc/bcell_from_all_datasets/obj.all.four_marker.csv")
obj.all.four_marker=(read.csv("E:/thesis/MSc/msc_figure/bcell_from_all_datasets/bcell_from_all_datasets/obj.all.four_marker.csv"))

View(obj.all.four_marker)

obj.all.four_marker = obj.all.four_marker %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

memory_gene=obj.all.four_marker[1:20,]
naive_gene=obj.all.four_marker[21:40,]
plasma_gene=obj.all.four_marker[41:60,]

memory_gene=obj.all.four_marker[1:50,]
naive_gene=obj.all.four_marker[51:100,]
plasma_gene=obj.all.four_marker[61:150,]

memory_gene=obj.all.four_marker[1:10,]
naive_gene=obj.all.four_marker[11:20,]
plasma_gene=obj.all.four_marker[21:30,]

memory_gene=obj.all.four_marker[1:30,]
naive_gene=obj.all.four_marker[11:20,]
plasma_gene=obj.all.four_marker[21:30,]
LIHC_mrna <- read_excel("E:/thesis/MSc/msc_figure/TCGA/cibersortx/LIHC_mrna.xlsx")


#coeffient

memoryGene_coeff=filter(LIHC_mrna,`Updated Name`%in% memory_gene$gene)
naiveGene_coeff=filter(LIHC_mrna,`Updated Name`%in% naive_gene$gene)
plasmaGene_coeff=filter(LIHC_mrna,`Updated Name`%in% plasma_gene$gene)

memoryGene_coeff_v=mean(memoryGene_coeff$`Cox coefficient`[1:30])
naiveGene_coeff_v=mean(naiveGene_coeff$`Cox coefficient`[1:30])
plasmaGene_coeff_v=mean(plasmaGene_coeff$`Cox coefficient`[1:30])

sum(memoryGene_coeff$`Cox coefficient`)
memoryCoeff=memoryGene_coeff[,c("Cox coefficient","Updated Name")]
memoryCoeff$subtype= "Memory Bcell"

naiveCoeff=naiveGene_coeff[,c("Cox coefficient","Updated Name")]
naiveCoeff$subtype= "Naive Bcell"

plasmaCoeff=plasmaGene_coeff[,c("Cox coefficient","Updated Name")]
plasmaCoeff$subtype= "Plasma cell"

df= rbind(memoryCoeff[1:30,],naiveCoeff[1:30,],plasmaCoeff[1:30,])
df$coefficient=df$`Cox coefficient`
library(ggpubr)



data(mpg)
ggplot(df, aes(x=subtype, y=coefficient, group=factor(subtype))) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(height = .3, width = .1))





ggbarplot(
  df, x = "subtype", y = "coefficient",
  add = c("mean_sd"),
  fill = c("#807F7F"),ylab="Cox coefficient",xlab="",width = 0.4
)+geom_jitter(position = position_jitter(height = .1, width = .08),size=1)+ylim(-0.4,0.4)+
  theme(axis.text.x = element_text(angle = 20,hjust = 0.5,vjust=0.4))+geom_hline(yintercept=0)




#checkpoint inhibitors gene

ics=c("TIGIT", "LY9", "PDCD1", 'LAG3', "CTLA4", "CD276", "NT5E", "PDCD1LG2", "CD274", "IDO1", "VSIR", "HAVCR2", "ENTPD", "CD274", "CTLA4", "HAVCR2", "LAG3", "TIGIT",
"PDCD1", "PDCD1LG2","SIGLEC15")
ics=unique(ics) #14
