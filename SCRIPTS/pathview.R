   ## LOAD LIBRARIES NEEDED

library(pathview)

## SET WORK DIRECTORY
setwd("/datos_2/FANCONI")


#>> PLOT FANCONI ANEMIA PATHWAY GENES IN EACH TYPE OF ANALYSIS 
celltype<-c("LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","MEP","Erythroid","Basophils") 



### mean values of all the data (fanconi and healthy)

for(i in 1:length(celltype)){
setwd("/datos_2/FANCONI/pathview/Mean_values")
  DE_genes_2008<-read.table(paste0("/datos_2/FANCONI/Fanconi_2008/Corrected_vs_uncorrected_",celltype[i],".txt"))
  logFC_2008<-as.data.frame(DE_genes_2008[,paste0("avg_logFC_" ,celltype[i])])
  rownames(logFC_2008)<-rownames(DE_genes_2008)
  colnames(logFC_2008) <- paste0("2008_",celltype[i])
  DE_genes_2006<-read.table(paste0("/datos_2/FANCONI/Fanconi_sample/DEA_all/Corrected_vs_uncorrected_",celltype[i],".txt"))
  logFC_2006<-as.data.frame(DE_genes_2006[,paste0("avg_logFC_" ,celltype[i])])
  rownames(logFC_2006)<-rownames(DE_genes_2006)
  colnames(logFC_2006) <- paste0("2006_",celltype[i])
  DE_genes_2004<-read.table(paste0("/datos_2/FANCONI/Fanconi_eltrombopag/Corrected_vs_uncorrected_",celltype[i],".txt"))
  logFC_2004<-as.data.frame(DE_genes_2004[,paste0("avg_logFC_" ,celltype[i])])
  rownames(logFC_2004)<-rownames(DE_genes_2004)
  colnames(logFC_2004) <- paste0("2004_",celltype[i])
  DE_genes_integration_2004<-read.table(paste0("/datos_2/FANCONI/Integration_eltrombopag/Healthy_vs_uncorrected_",celltype[i],".txt"))
  logFC_integration_2004<-as.data.frame(DE_genes_integration_2004[,paste0("avg_logFC_" ,celltype[i])])
  rownames(logFC_integration_2004)<-rownames(DE_genes_integration_2004)
  colnames(logFC_integration_2004) <- paste0("2004_integration_",celltype[i])
  
    DE_genes_integration_2006<-read.table(paste0("/datos_2/FANCONI/Integration_new_healthy/DEA_all/Healthy_vs_uncorrected_",celltype[i],".txt"))
  logFC_integration_2006<-as.data.frame(DE_genes_integration_2006[,paste0("avg_logFC_" ,celltype[i])])
  rownames(logFC_integration_2006)<-rownames(DE_genes_integration_2006)
  colnames(logFC_integration_2006) <- paste0("2006_integration_",celltype[i])
  
  DE_genes_integration_2008<-read.table(paste0("/datos_2/FANCONI/Fanconi_2008/integration/Healthy_uncorrected_",celltype[i],".txt"))
  logFC_integration_2008<-as.data.frame(DE_genes_integration_2008[,paste0("avg_logFC_" ,celltype[i])])
  rownames(logFC_integration_2008)<-rownames(DE_genes_integration_2008)
  colnames(logFC_integration_2008) <- paste0("2008_integration_",celltype[i])
  
MyMerge       <- function(x, y){
  df            <- merge(x, y, by= "row.names", all.x= T, all.y= T)
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)
}
dat           <- Reduce(MyMerge, list(logFC_2004, logFC_2006, logFC_2008, logFC_integration_2004,logFC_integration_2006,logFC_integration_2008))
dat[is.na(dat)]<-0

#dat$fanconi<-rowMeans(dat[,1:3])
#dat$integration<-rowMeans(dat[,4:6])
#
#dat<-dat[,7:8]
dat[is.na(dat)]<-0
  #png(paste0("Fanconi_sample/pathview/DNA_replication_pathway_",celltype[i],".png"))
  p1<-pathview(gene.data = dat, species = "hsa", pathway.id = "05014", gene.idtype = "SYMBOL", out.suffix =celltype[i])
    #dev.off()
  }



## only healthy information (healthy vs uncorrected)

for(i in 1:length(celltype)){
setwd("/datos_2/FANCONI/pathview/Healthy_vs_uncorrected")
  
  DE_genes_integration_2004<-read.table(paste0("/datos_2/FANCONI/Integration_eltrombopag/Healthy_vs_uncorrected_",celltype[i],".txt"))
  logFC_integration_2004<-as.data.frame(DE_genes_integration_2004[,paste0("avg_logFC_" ,celltype[i])])
  rownames(logFC_integration_2004)<-rownames(DE_genes_integration_2004)
  colnames(logFC_integration_2004) <- paste0("2004_integration_",celltype[i])
  
    DE_genes_integration_2006<-read.table(paste0("/datos_2/FANCONI/Integration_new_healthy/DEA_all/Healthy_vs_uncorrected_",celltype[i],".txt"))
  logFC_integration_2006<-as.data.frame(DE_genes_integration_2006[,paste0("avg_logFC_" ,celltype[i])])
  rownames(logFC_integration_2006)<-rownames(DE_genes_integration_2006)
  colnames(logFC_integration_2006) <- paste0("2006_integration_",celltype[i])
  
  DE_genes_integration_2008<-read.table(paste0("/datos_2/FANCONI/Fanconi_2008/integration/Healthy_uncorrected_",celltype[i],".txt"))
  logFC_integration_2008<-as.data.frame(DE_genes_integration_2008[,paste0("avg_logFC_" ,celltype[i])])
  rownames(logFC_integration_2008)<-rownames(DE_genes_integration_2008)
  colnames(logFC_integration_2008) <- paste0("2008_integration_",celltype[i])
  
MyMerge       <- function(x, y){
  df            <- merge(x, y, by= "row.names", all.x= T, all.y= T)
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)
}
dat           <- Reduce(MyMerge, list(logFC_integration_2004,logFC_integration_2006,logFC_integration_2008))
dat[is.na(dat)]<-0

  #png(paste0("Fanconi_sample/pathview/DNA_replication_pathway_",celltype[i],".png"))
  p1<-pathview(gene.data = dat, species = "hsa", pathway.id = "03460", gene.idtype = "SYMBOL", out.suffix =celltype[i], same.layer=F)
    #dev.off()
  }
  
  
## only fanconi information (healthy vs uncorrected)

for(i in 1:length(celltype)){
setwd("/datos_2/FANCONI/pathview")

  DE_genes_2008<-read.table(paste0("/datos_2/FANCONI/Fanconi_2008/Corrected_vs_uncorrected_",celltype[i],".txt"))
  logFC_2008<-as.data.frame(DE_genes_2008[,paste0("avg_logFC_" ,celltype[i])])
  rownames(logFC_2008)<-rownames(DE_genes_2008)
  colnames(logFC_2008) <- paste0("2008_",celltype[i])
  DE_genes_2006<-read.table(paste0("/datos_2/FANCONI/Fanconi_sample/DEA_all/Corrected_vs_uncorrected_",celltype[i],".txt"))
  logFC_2006<-as.data.frame(DE_genes_2006[,paste0("avg_logFC_" ,celltype[i])])
  rownames(logFC_2006)<-rownames(DE_genes_2006)
  colnames(logFC_2006) <- paste0("2006_",celltype[i])
  DE_genes_2004<-read.table(paste0("/datos_2/FANCONI/Fanconi_eltrombopag/Corrected_vs_uncorrected_",celltype[i],".txt"))
  logFC_2004<-as.data.frame(DE_genes_2004[,paste0("avg_logFC_" ,celltype[i])])
  rownames(logFC_2004)<-rownames(DE_genes_2004)
  colnames(logFC_2004) <- paste0("2004_",celltype[i])
  
  
MyMerge       <- function(x, y){
  df            <- merge(x, y, by= "row.names", all.x= T, all.y= T)
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)
}
dat           <- Reduce(MyMerge, list(logFC_2004, logFC_2006, logFC_2008))
dat[is.na(dat)]<-0

  #png(paste0("Fanconi_sample/pathview/DNA_replication_pathway_",celltype[i],".png"))
  p1<-pathview(gene.data = dat, species = "hsa", pathway.id = "04930", gene.idtype = "SYMBOL", out.suffix =celltype[i])
    #dev.off()
  }