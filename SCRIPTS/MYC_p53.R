library(Seurat)
library(biomaRt)
library(ggplot2)
library(clusterProfiler)
library('org.Hs.eg.db') 
library(ggpubr)
library(patchwork) 
library(pathview)
library(xlsx)
set.seed(1234567)
setwd("/datos_2/FANCONI")

## read pathway information


celltype<-c("HSC","LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","MEP","Erythroid","Basophils")
gmtfile <- read.gmt('Fanconi_02006/MYC_p53/MYC_p53_geneset.gmt.txt')

for (i in 1:length(celltype)){
  set.seed(1234567)
  dea_res<-read.table(paste0("Integration_eltrombopag/Healthy_vs_uncorrected_",celltype[i],".txt"))
  
  
  ranked_genes<-dea_res[,paste0("avg_logFC_" ,celltype[i])]
  names(ranked_genes)<-rownames(dea_res)
  ranked_genes <- sort(ranked_genes, decreasing = T)
  head(ranked_genes)
  tail(ranked_genes)
  
  egmt_Basophils_uncorrected_2004<- GSEA(ranked_genes, TERM2GENE=gmtfile, verbose=T, minGSSize = 10, pvalueCutoff = 1)
  pdf(paste0("Integration_eltrombopag/MYC_p53/ridgeplot_",celltype[i],".pdf"))
  ridgeplot(egmt_Basophils_uncorrected_2004)
  dev.off()
  
  
 
 }
 