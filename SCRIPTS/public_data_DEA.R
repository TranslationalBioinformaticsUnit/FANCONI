




library(Seurat)
library(ggplot2)
library(clusterProfiler)


setwd("/datos_2/FANCONI/Public_data")

integration<-readRDS("integrated_cellnames.RDS")

integration[["percent.mt"]] <- PercentageFeatureSet(integration, pattern = "^MT-")
integration <- subset(integration, subset = nFeature_RNA > 2000 & nCount_RNA < 8e04 & percent.mt < 4)


integration$CellType<-integration$predicted.id
integration@meta.data$CellType<-factor(x=integration@meta.data$CellType, levels=c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils"))

integration<-SetIdent(integration, value="CellType")

levels(x = integration)<-c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils")


integration<-NormalizeData(integration)
integration <- ScaleData(integration, features = rownames(integration))


integration@meta.data$CD34 <- integration@assays$RNA@data["CD34",]



## filter CD34 positive cells
integration<-subset(integration, subset=CD34>0.125)

integration@meta.data$FANCA <- integration@assays$RNA@data["FANCA",]
integration@meta.data$MYC <- integration@assays$RNA@data["MYC",]


# create a new object with myc gene expression filtered
MYC_filtered<-subset(integration, subset=MYC>0)

# plot MYC expression by cell type, and fill by project
pdf("Analysis/MYC_plots/MYC_expression_cell_type.pdf", width = 10)
ggplot(fanca_filtered@meta.data, aes(x=Sample, y=FANCA, fill=Sample)) + 
  geom_boxplot()+
  labs(title = "MYC expression by Sample")

dev.off()



fanca_filtered<-subset(integration, subset=FANCA>0)
MYC_filtered<-subset(integration, subset=MYC>0)

integration@meta.data$MYC <- integration@assays$RNA@data["MYC",]


##Clustering

DefaultAssay(integration)<-"integrated"

integration <- RunPCA(integration, npcs = 50, verbose = FALSE)
P2<-DimPlot(integration)
P2

# Select components to use in further steps
ElbowPlot(integration,ndims = 50)

## UMAP
integration <- RunUMAP(integration, reduction = "pca", dims = 1:8)

## Clustering
integration <- FindNeighbors(integration, reduction = "pca", dims = 1:8)
integration <- FindClusters(integration, resolution = 0.8, n.start = 1000) ## n.start=1000 para obtener resultados similares dependiendo del orden de células

# Transform Numbers to Letters
integration$Cluster<-LETTERS[as.numeric(as.vector(integration$seurat_clusters))+1]
LX<-length(levels(integration$seurat_clusters))
integration$Cluster<-factor(integration$Cluster,levels = LETTERS[1:LX])
integration<-SetIdent(integration,value = "Cluster")

## DimPlot
pdf("Integration/umap_clustering.pdf")
DimPlot(integration,label=TRUE,label.size = 7)
dev.off()

##label transfer
healthy<-readRDS("/datos_2/FANCONI/DATA/seurat_young.rds")

DefaultAssay(integration)<-"integrated"
DEfaultAssay(healthy)<-"integrated"

fanconi.anchors <- FindTransferAnchors(reference = healthy, query = integration, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = fanconi.anchors, refdata = healthy$CellType, 
                            dims = 1:30)
integration <- AddMetaData(integration, metadata = predictions)


#DEA

Cell_type<-c("LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid", "Basophils") #HSC cab't be included, any corrected cell

for(i in 1:length(Cell_type)) {
  
disease_markers<-FindMarkers(integration, ident.1="Healthy", ident.2="Fanconi", group.by= "Project", logfc.threshold = 0, subset.ident = Cell_type[i])
sig_genes_up <- which(disease_markers$p_val_adj<0.05 & disease_markers$avg_logFC>0.25)
sig_genes_down <- which(disease_markers$p_val_adj<0.05 & disease_markers$avg_logFC< -0.25)
disease_markers$sig_genes <- rep(0,nrow(disease_markers))
  disease_markers$sig_genes[disease_markers$avg_logFC>0] <-0.5
  disease_markers$sig_genes[disease_markers$avg_logFC<0] <- -0.5
  disease_markers$sig_genes[sig_genes_up] <- 1
  disease_markers$sig_genes[sig_genes_down] <- -1
  colnames(disease_markers)<-paste(colnames(disease_markers), Cell_type[i], sep = "_")
write.table(disease_markers, paste0("Analysis/DEA/Healthy_fanconi_",Cell_type[i],".txt", sep=""),sep="\t",quote=F,row.names=T)
markers<-NULL
}



# GSEA
Cell_type<-c("LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","MEP","Erythroid", "Basophils")
gmtfile <- read.gmt('/datos_2/FANCONI/Fanconi_sample/GSEA_public_paper/MYC_pathway_genes.gmt.txt')


for (i in 1:length(Cell_type)){
  set.seed(1234567)
  dea_res<-read.table(paste0("Healthy_uncorrected_",Cell_type[i],".txt"))
  
  
  ranked_genes<-dea_res[,paste0("avg_logFC_" ,Cell_type[i])]
  names(ranked_genes)<-rownames(dea_res)
  ranked_genes <- sort(ranked_genes, decreasing = T)
  head(ranked_genes)
  tail(ranked_genes)
  
  egmt<- GSEA(ranked_genes, TERM2GENE=gmtfile, verbose=T, minGSSize = 10, pvalueCutoff = 1)
  pdf(paste0("MYC_pathway/ridgeplot_",Cell_type[i],".pdf"))
  ridgeplot(egmt)
  dev.off()
  
  
 
 }

