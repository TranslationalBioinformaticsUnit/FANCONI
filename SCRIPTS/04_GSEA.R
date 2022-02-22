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

pathways<-read.xlsx("Pathways_of_interest.xlsx", sheetName="Hoja1")

pathways<-as.character(na.omit(pathways$ridgeplot))



celltype<-c("HSC","LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils")


for (i in 1:length(celltype)){
  set.seed(1234567)
  dea_res<-read.table(paste0("Fanconi_2002/Integration/Healthy_Uncorrected_",celltype[i],".txt")) # load the results of DEA analysis
  
  ranked_genes<-dea_res[,paste0("avg_logFC_" ,celltype[i])] # take only the value of the logFC
  names(ranked_genes)<-rownames(dea_res) # order the genes depending on the logFC value
    ranked_genes <- sort(ranked_genes, decreasing = T)
  head(ranked_genes)
  tail(ranked_genes)
  
  
  
  gse_Basohphils_integration_2002<-gseGO(geneList = ranked_genes, OrgDb = org.Hs.eg.db, #calcula the enrichment with GO pathways
             keyType = 'SYMBOL',
             ont = "BP", nPerm = 1000,
             minGSSize = 10, maxGSSize = 200,
             pvalueCutoff = 1,
             verbose=FALSE)
             
  gse_LMPP_uncorrected_2006_2<-gse_LMPP_uncorrected_2006
  prueba<-gse_LMPP_uncorrected_2006[gse_LMPP_uncorrected_2006$Description %in% pathways,]
  gse_Basophils_2008_2@result<-prueba
  pdf("Integration_new_healthy/NEW_GSEA/LMPP_uncorrected_ridgeplot_filtered.pdf", width=15)
  ridgeplot(gse_Basophils_2008_2, showCategory =length(pathways))
  dev.off()

  write.table(gse_Basohphils_integration_2002,paste0("Fanconi_2002/Integration/GSEA/GO_",celltype[i],".txt", sep=""),
  sep="\t",quote=F,row.names=F) # write the results of the analysis 
 
}

## Ridgeplot of all the data together, to understand the function check our_ridgeplot.R script

source("/datos_2/FANCONI/SCRIPTS/our_ridgeplot.R")

#How to run:
#list_gsea <- list(gse_LMPP,gse_Cycling_LMPP, gse_GMP1, gse_GMP2, gse_Monocytes, gse_DC, gse_CLP,gse_PreB,gse_MEP, gse_Erythroid, gse_Basophils)

list_gsea <- list(gse_Monocytes_2002, gse_Monocytes_integration_2002, gse_Monocytes_2004, gse_Monocytes_uncorrected_2004, gse_Monocytes, gse_Monocytes_uncorrected_2006, gse_Monocytes_2008,gse_Monocytes_uncorrected_2008)
#gsea_names <- c("HSC","LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid", "Basophils")

gsea_names <- c("Corrected_uncorrected_2002","Healthy_uncorrected_2002","Corrected_uncorrected_2004","Healthy_uncorrected_2004","Corrected_uncorrected_2006","Healthy_uncorrected_2006","Corrected_uncorrected_2008", "Healthy_uncorrected_2008")
#order_by <- "Corrected_uncorrected_2008"

my_output_2 <- our_ridge_multiplot(list_gsea, path, gsea_names, order_by=NULL)





