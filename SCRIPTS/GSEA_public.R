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
setwd("/datos_2/FANCONI/Public_data")

## read pathway information

celltype<-c("HSC","LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","MEP","Erythroid")


for (i in 1:length(celltype)){
  set.seed(1234567)
  dea_res<-read.table(paste0("Analysis/DEA/Healthy_fanconi_",celltype[i],".txt"))
  
  
  ranked_genes<-dea_res[,paste0("avg_logFC_" ,celltype[i])]
  names(ranked_genes)<-rownames(dea_res)
  ranked_genes <- sort(ranked_genes, decreasing = T)
  head(ranked_genes)
  tail(ranked_genes)
  
  
  
  gse_Erythroid_public<-gseGO(geneList = ranked_genes, OrgDb = org.Hs.eg.db,
             keyType = 'SYMBOL',
             ont = "BP", nPerm = 1000,
             minGSSize = 10, maxGSSize = 200,
             pvalueCutoff = 1,
             verbose=FALSE)
             
  write.table(gse_Erythroid_public,paste0("Analysis/GSEA/GO_uncorrected_",celltype[i],".txt", sep=""),
  sep="\t",quote=F,row.names=F)
 
}



source("/datos_2/FANCONI/SCRIPTS/our_ridgeplot.R")

#How to run:
list_gsea <- list(gse_HSC_public, gse_LMPP_public,gse_Cycling_LMPP_public, gse_GMP1_public, gse_GMP2_public, gse_Monocytes_public, gse_DC_public, gse_CLP_public, gse_MEP_public, gse_Erythroid_public)

#list_gsea <- list(gse_HSC_2004, gse_HSC_uncorrected_2004)
gsea_names <- c("HSC","LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","MEP","Erythroid")
#order_by <- "Corrected_uncorrected_2008"

my_output_2 <- our_ridge_multiplot(list_gsea, path, gsea_names, order_by=NULL)

