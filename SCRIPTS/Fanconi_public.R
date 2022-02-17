## Analysis of the public data

library(Seurat)
set.seed(1234567)

setwd("/datos_2/FANCONI/Public_data")

#load DATA
prueba<-readRDS("integration_berria.RDS")

# deffault assay RNA

DefaultAssay(prueba)<-"RNA"


prueba@meta.data$CD34<- prueba@assays$RNA@data["CD34",]

# Initial QC Netrics
# We get the % of UMIs mapped to mitochondrial genes as a way to identify "broken/stressed/ cells
## Reference
#Ilicic, T., Kim, J.K., Kolodziejczyk, A.A. et al. Classification of low quality cells from single-cell RNA-seq data. Genome Biol 17, 29 (2016). 
#https://doi.org/10.1186/s13059-016-0888-1
##
# Also the # of UMIs per cell and the # of Genes per cell are used
prueba[["percent.mt"]] <- PercentageFeatureSet(prueba, pattern = "^MT-")
VlnPlot(prueba,c("nCount_RNA","nFeature_RNA","percent.mt"))

# Filtering based on QC parameters
prueba <- subset(prueba, subset = nFeature_RNA > 2000  & nCount_RNA < 8e04 & percent.mt < 4)

# Normalize and Scale Data
prueba<-NormalizeData(prueba)
prueba <- ScaleData(prueba, features = rownames(prueba))

## Cell Cycle Analysis
prueba<-CellCycleScoring(prueba,s.features = cc.genes$s.genes,g2m.features = cc.genes$g2m.genes)

# Using PCA, we see if there is an effect caussed due to cell cycle phase
prueba<-RunPCA(prueba,features = c(cc.genes$s.genes,cc.genes$g2m.genes))
DimPlot(prueba,group.by = "Phase")
prueba$CC.Difference <- prueba$S.Score - prueba$G2M.Score

# Do the clustering in INTEGRATION slot


DefaultAssay(integration)<-"integrated"

## Clustering

# PCA
prueba <- RunPCA(prueba, npcs = 50, verbose = FALSE)
P1<-DimPlot(prueba,group.by = "Phase")
P2<-DimPlot(prueba)
P1+P2
# Select components to use in further steps
ElbowPlot(prueba,ndims = 50)
## UMAP
prueba <- RunUMAP(prueba, reduction = "pca", dims = 1:6)

## Clustering
prueba <- FindNeighbors(prueba, reduction = "pca", dims = 1:6)
prueba <- FindClusters(prueba, resolution = 0.8, n.start = 1000) ## n.start=1000 para obtener resultados similares dependiendo del orden de células

# Transform Numbers to Letters
integration$Cluster<-LETTERS[as.numeric(as.vector(integration$seurat_clusters))+1]
LX<-length(levels(integration$seurat_clusters))
integration$Cluster<-factor(integration$Cluster,levels = LETTERS[1:LX])
integration<-SetIdent(integration,value = "Cluster")

# Distribution of FANCA gene

integration_filt@meta.data$FANCA <- integration_filt@assays$RNA@data["FANCA",]

integration_FANCA_filtered<-subset(integration_filt, subset=FANCA>0)

pdf("RESULTS/FANCA_distr_violin_filtered.pdf")
VlnPlot(integration,"MYC")
dev.off()

ggplot(integration_FANCA_filtered@meta.data, aes(x=Sample, y=FANCA, fill=Sample)) + 
  geom_boxplot()+
  labs(title = "FANCA expression by Sample")

#Label transfer to name the cells

healthy<-readRDS("/datos_2/FANCONI/DATA/seurat_young.rds")

fanconi.anchors <- FindTransferAnchors(reference = healthy, query = prueba, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = fanconi.anchors, refdata = healthy$CellType, 
                            dims = 1:30)
prueba <- AddMetaData(prueba, metadata = predictions)


p1<- ggplot(prueba@meta.data, aes(x=predicted.id, fill=Sample)) + 
  geom_bar(position = "fill")+
  labs(title="% Cell type by sample",x = "Sample", y = "%")
ggsave(plot = p1, height = 7, width = 9, dpi=600, filename = "percentage_by_sample.pdf")


Cell_type<-c("HSC","LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils") #HSC cab't be included, any corrected cell

for(i in 1:length(Cell_type)) {
  
disease_markers<-FindMarkers(prueba, ident.1="Fanconi", ident.2="Healthy", group.by= "Project", logfc.threshold = 0.25, subset.ident = Cell_type[i])
write.table(disease_markers, paste0("Fanconi_vs_healthy_",Cell_type[i],".txt", sep=""),sep="\t",quote=F,row.names=T)
markers<-NULL
}

