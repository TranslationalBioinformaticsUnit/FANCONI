

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