#################### Genes obtained from DEA by cell type ####################

set.seed(1234567)
##set work directory
setwd("datos_2/FANCONI")

library(GeneSetCluster)   # https://github.com/TranslationalBioinformaticsUnit/GeneSetCluster/wiki




# ---------------------

# STEP 1.1: LOAD THE DATA 

##Fanconi patient

GO_LMPP_corrected_all<-read.table("Integration_updated/GSEA_all_corrected/GO/GO_LMPP.txt", sep="\t", header=T)
GO_Cycling_LMPP_corrected_all<-read.table("Integration_updated/GSEA_all_corrected/GO/GO_Cycling_LMPP.txt", sep="\t", header=T)
GO_GMP1_corrected_all<-read.table("Integration_updated/GSEA_all_corrected/GO/GO_Cycling_LMPP.txt", sep="\t", header=T)
GO_GMP2_corrected_all<-read.table("Integration_updated/GSEA_all_corrected/GO/GO_Cycling_LMPP.txt", sep="\t", header=T)
GO_Monocytes_corrected_all<-read.table("Integration_updated/GSEA_all_corrected/GO/GO_Cycling_LMPP.txt", sep="\t", header=T)
GO_DC_corrected_all<-read.table("Integration_updated/GSEA_all_corrected/GO/GO_Cycling_LMPP.txt", sep="\t", header=T)
GO_CLP_corrected_all<-read.table("Integration_updated/GSEA_all_corrected/GO/GO_Cycling_LMPP.txt", sep="\t", header=T)
GO_PreB_corrected_all<-read.table("Integration_updated/GSEA_all_corrected/GO/GO_Cycling_LMPP.txt", sep="\t", header=T)
GO_MEP_corrected_all<-read.table("Integration_updated/GSEA_all_corrected/GO/GO_Cycling_LMPP.txt", sep="\t", header=T)
GO_Erythroid_corrected_all<-read.table("Integration_updated/GSEA_all_corrected/GO/GO_Cycling_LMPP.txt", sep="\t", header=T)
GO_Basophils_corrected_all<-read.table("Integration_updated/GSEA_all_corrected/GO/GO_Cycling_LMPP.txt", sep="\t", header=T)


## take only the significant genes in each case

dea_LMPP_fanconi<-rownames(dea_LMPP_fanconi[dea_LMPP_fanconi$p_val_adj<0.05,])
dea_Cycling_LMPP_fanconi<-rownames(dea_Cycling_LMPP_fanconi[dea_Cycling_LMPP_fanconi$p_val_adj<0.05,])
dea_GMP1_fanconi<-rownames(dea_GMP1_fanconi[dea_GMP1_fanconi$p_val_adj<0.05,])
dea_GMP2_fanconi<-rownames(dea_GMP2_fanconi[dea_GMP2_fanconi$p_val_adj<0.05,])
dea_Monocytes_fanconi<-rownames(dea_Monocytes_fanconi[dea_Monocytes_fanconi$p_val_adj<0.05,])
dea_DC_fanconi<-rownames(dea_DC_fanconi[dea_DC_fanconi$p_val_adj<0.05,])
dea_CLP_fanconi<-rownames(dea_CLP_fanconi[dea_CLP_fanconi$p_val_adj<0.05,])
dea_PreB_fanconi<-rownames(dea_PreB_fanconi[dea_PreB_fanconi$p_val_adj<0.05,])
dea_MEP_fanconi<-rownames(dea_MEP_fanconi[dea_MEP_fanconi$p_val_adj<0.05,])
dea_Erythroid_fanconi<-rownames(dea_Erythroid_fanconi[dea_Erythroid_fanconi$p_val_adj<0.05,])
dea_Basophils_fanconi<-rownames(dea_Basophils_fanconi[dea_Basophils_fanconi$p_val_adj<0.05,])


## take all the genes

dea_LMPP_fanconi<-rownames(dea_LMPP_fanconi)
dea_Cycling_LMPP_fanconi<-rownames(dea_Cycling_LMPP_fanconi)
dea_GMP1_fanconi<-rownames(dea_GMP1_fanconi)
dea_GMP2_fanconi<-rownames(dea_GMP2_fanconi)
dea_Monocytes_fanconi<-rownames(dea_Monocytes_fanconi)
dea_DC_fanconi<-rownames(dea_DC_fanconi)
dea_CLP_fanconi<-rownames(dea_CLP_fanconi)
dea_PreB_fanconi<-rownames(dea_PreB_fanconi)
dea_MEP_fanconi<-rownames(dea_MEP_fanconi)
dea_Erythroid_fanconi<-rownames(dea_Erythroid_fanconi)
dea_Basophils_fanconi<-rownames(dea_Basophils_fanconi)


## Healthy donor

dea_HSC_healthy<-read.table("TABLES/Healthy/markers_fanconi_positive_vs_negative_HSC.txt", sep="\t", header=T)
dea_LMPP_healthy<-read.table("TABLES/Healthy/markers_fanconi_positive_vs_negative_LMPP.txt", sep="\t", header=T)
dea_Cycling_LMPP_healthy<-read.table("TABLES/Healthy/markers_fanconi_positive_vs_negative_Cycling_LMPP.txt", sep="\t", header=T)
dea_GMP1_healthy<-read.table("TABLES/Healthy/markers_fanconi_positive_vs_negative_GMP1.txt", sep="\t", header=T)
dea_GMP2_healthy<-read.table("TABLES/Healthy/markers_fanconi_positive_vs_negative_GMP2.txt", sep="\t", header=T)
dea_Monocytes_healthy<-read.table("TABLES/Healthy/markers_fanconi_positive_vs_negative_Monocytes.txt", sep="\t", header=T)
dea_DC_healthy<-read.table("TABLES/Healthy/markers_fanconi_positive_vs_negative_DC.txt", sep="\t", header=T)
dea_CLP_healthy<-read.table("TABLES/Healthy/markers_fanconi_positive_vs_negative_CLP.txt", sep="\t", header=T)
dea_PreB_healthy<-read.table("TABLES/Healthy/markers_fanconi_positive_vs_negative_PreB.txt", sep="\t", header=T)
dea_MEP_healthy<-read.table("TABLES/Healthy/markers_fanconi_positive_vs_negative_MEP.txt", sep="\t", header=T)
dea_Erythroid_healthy<-read.table("TABLES/Healthy/markers_fanconi_positive_vs_negative_Erythroid.txt", sep="\t", header=T)
dea_Basophils_healthy<-read.table("TABLES/Healthy/markers_fanconi_positive_vs_negative_Basophils.txt", sep="\t", header=T)

## take the healthy significants

dea_HSC_healthy<-rownames(dea_HSC_healthy[dea_HSC_healthy$p_val_adj<0.05,])
dea_LMPP_healthy<-rownames(dea_LMPP_healthy[dea_LMPP_healthy$p_val_adj<0.05,])
dea_Cycling_LMPP_healthy<-rownames(dea_Cycling_LMPP_healthy[dea_Cycling_LMPP_healthy$p_val_adj<0.05,])
dea_GMP1_healthy<-rownames(dea_GMP1_healthy[dea_GMP1_healthy$p_val_adj<0.05,])
dea_GMP2_healthy<-rownames(dea_GMP2_healthy[dea_GMP2_healthy$p_val_adj<0.05,])
dea_Monocytes_healthy<-rownames(dea_Monocytes_healthy[dea_Monocytes_healthy$p_val_adj<0.05,])
dea_DC_healthy<-rownames(dea_DC_healthy[dea_DC_healthy$p_val_adj<0.05,])
dea_CLP_healthy<-rownames(dea_CLP_healthy[dea_CLP_healthy$p_val_adj<0.05,])
dea_PreB_healthy<-rownames(dea_PreB_healthy[dea_PreB_healthy$p_val_adj<0.05,])
dea_MEP_healthy<-rownames(dea_MEP_healthy[dea_MEP_healthy$p_val_adj<0.05,])
dea_Erythroid_healthy<-rownames(dea_Erythroid_healthy[dea_Erythroid_healthy$p_val_adj<0.05,])
dea_Basophils_healthy<-rownames(dea_Basophils_healthy[dea_Basophils_healthy$p_val_adj<0.05,])

## take all the DE genes 

dea_HSC_healthy<-rownames(dea_HSC_healthy)
dea_LMPP_healthy<-rownames(dea_LMPP_healthy)
dea_Cycling_LMPP_healthy<-rownames(dea_Cycling_LMPP_healthy)
dea_GMP1_healthy<-rownames(dea_GMP1_healthy)
dea_GMP2_healthy<-rownames(dea_GMP2_healthy)
dea_Monocytes_healthy<-rownames(dea_Monocytes_healthy)
dea_DC_healthy<-rownames(dea_DC_healthy)
dea_CLP_healthy<-rownames(dea_CLP_healthy)
dea_PreB_healthy<-rownames(dea_PreB_healthy)
dea_MEP_healthy<-rownames(dea_MEP_healthy)
dea_Erythroid_healthy<-rownames(dea_Erythroid_healthy)
dea_Basophils_healthy<-rownames(dea_Basophils_healthy)


MyMerge       <- function(x, y){
  df            <- merge(x, y, by= "row.names", all.x= T, all.y= T)
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)
}


dea_all <- Reduce(MyMerge, list(dea_LMPP_fanconi, dea_Cycling_LMPP_fanconi, dea_GMP1_fanconi, dea_GMP2_fanconi,dea_Monocytes_fanconi,dea_DC_fanconi,dea_CLP_fanconi,dea_PreB_fanconi,dea_MEP_fanconi,dea_Erythroid_fanconi,dea_Basophils_fanconi,dea_HSC_healthy,dea_LMPP_healthy,dea_Cycling_LMPP_healthy,dea_GMP1_healthy,dea_GMP2_healthy,dea_Monocytes_healthy,dea_DC_healthy,dea_DC_healthy,dea_CLP_healthy,dea_PreB_healthy,dea_MEP_healthy,dea_Erythroid_healthy,dea_Basophils_healthy))

background <- unique(c(dea_LMPP_fanconi, dea_Cycling_LMPP_fanconi, dea_GMP1_fanconi, dea_GMP2_fanconi,dea_Monocytes_fanconi,dea_DC_fanconi,dea_CLP_fanconi,dea_PreB_fanconi,dea_MEP_fanconi,dea_Erythroid_fanconi,dea_Basophils_fanconi,dea_HSC_healthy,dea_LMPP_healthy,dea_Cycling_LMPP_healthy,dea_GMP1_healthy,dea_GMP2_healthy,dea_Monocytes_healthy,dea_DC_healthy,dea_DC_healthy,dea_CLP_healthy,dea_PreB_healthy,dea_MEP_healthy,dea_Erythroid_healthy,dea_Basophils_healthy))

#dea_all <- rbind(dea_LMPP_fanconi,dea_Cycling_LMPP_fanconi,dea_GMP1_fanconi,dea_GMP2_fanconi,dea_Monocytes_fanconi,dea_DC_fanconi,dea_CLP_fanconi,dea_PreB_fanconi,dea_MEP_fanconi,dea_Erythroid_fanconi,dea_Basophils_fanconi,
#dea_HSC_healthy,dea_LMPP_healthy,dea_Cycling_LMPP_healthy,dea_GMP1_healthy,dea_GMP2_healthy,dea_Monocytes_healthy,dea_DC_healthy,dea_CLP_healthy,dea_PreB_healthy,dea_MEP_healthy,dea_Erythroid_healthy,dea_Basophils_healthy)

fanconi_genes<-c("BRCA2", "BRIP1", "ERCC4", "FANCA", "FANCB", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI", "FANCL", "FANCM", "PALB2", "RAD51C", "SLX4", "UBE2T", "XRCC2", "XPF", "RD51", "REV7","RFWD3")


# STEP 1.2: MERGE THE DATA

dea <- as.data.frame(background)

rownames(dea)<-dea$background
# Fanconi sample

dea$LMPP_fanconi <- ifelse(rownames(dea) %in% dea_LMPP_fanconi, "yes", "no")
dea$cycling_LMPP_fanconi <- ifelse(rownames(dea) %in% dea_Cycling_LMPP_fanconi, "yes", "no")
dea$GMP1_fanconi <- ifelse(rownames(dea) %in% dea_GMP1_fanconi, "yes", "no")
dea$GMP2_fanconi <- ifelse(rownames(dea) %in% dea_GMP2_fanconi, "yes", "no")
dea$Monocytes_fanconi <- ifelse(rownames(dea) %in% dea_Monocytes_fanconi, "yes", "no")
dea$DC_fanconi <- ifelse(rownames(dea) %in% dea_DC_fanconi, "yes", "no")
dea$CLP_fanconi<- ifelse(rownames(dea) %in% dea_CLP_fanconi, "yes", "no")
dea$PreB_fanconi <- ifelse(rownames(dea) %in% dea_PreB_fanconi, "yes", "no")
dea$MEP_fanconi <- ifelse(rownames(dea) %in% dea_MEP_fanconi, "yes", "no")
dea$erythroid_fanconi <- ifelse(rownames(dea) %in% dea_Erythroid_fanconi, "yes", "no")
dea$basophils_fanconi<- ifelse(rownames(dea) %in% dea_Basophils_fanconi, "yes", "no")

#Healthy donor

dea$HSC_healthy <- ifelse(rownames(dea) %in% dea_HSC_healthy, "yes", "no")
dea$LMPP_healthy <- ifelse(rownames(dea) %in% dea_LMPP_healthy, "yes", "no")
dea$cycling_LMPP_healthy <- ifelse(rownames(dea) %in% dea_Cycling_LMPP_healthy, "yes", "no")
dea$GMP1_healthy <- ifelse(rownames(dea) %in% dea_GMP1_healthy, "yes", "no")
dea$GMP2_healthy <- ifelse(rownames(dea) %in% dea_GMP2_healthy, "yes", "no")
dea$Monocytes_healthy <- ifelse(rownames(dea) %in% dea_Monocytes_healthy, "yes", "no")
dea$DC_healthy <- ifelse(rownames(dea) %in% dea_DC_healthy, "yes", "no")
dea$CLP_healthy<- ifelse(rownames(dea) %in% dea_CLP_healthy, "yes", "no")
dea$PreB_healthy <- ifelse(rownames(dea) %in% dea_PreB_healthy, "yes", "no")
dea$MEP_healthy <- ifelse(rownames(dea) %in% dea_MEP_healthy, "yes", "no")
dea$erythroid_healthy <- ifelse(rownames(dea) %in% dea_Erythroid_healthy, "yes", "no")
dea$basophils_healthy<- ifelse(rownames(dea) %in% dea_Basophils_healthy, "yes", "no")
dea$fanconi_gene<-ifelse(rownames(dea)%in% fanconi_genes, "yes","no")



head(dea)