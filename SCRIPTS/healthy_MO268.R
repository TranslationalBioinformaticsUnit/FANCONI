###### analysis of the healthy 268 sample


# load required libraries
library(Seurat)
library(ggplot2)
library(ggpubr)

set.seed(1234567)
##set work directory
setwd("/datos_2/FANCONI")
#####
# CD34 Bone Marrow

# Read CellRanger Output (.h5 files are faster to load)
mo268<-Read10X("DATA/filtered_feature_bc_matrix")

# Create Seurat Object with initial filtering
# A gene must be expressed in at  least 3 cells to be considered
# A cell must have at least 400 detected genes to be considered
mo268 <- CreateSeuratObject(counts = mo268, project = "Healthy", min.cells = 3, min.features = 400)

# Initial QC Netrics
# We get the % of UMIs mapped to mitochondrial genes as a way to identify "broken/stressed/ cells
## Reference
#Ilicic, T., Kim, J.K., Kolodziejczyk, A.A. et al. Classification of low quality cells from single-cell RNA-seq data. Genome Biol 17, 29 (2016). 
#https://doi.org/10.1186/s13059-016-0888-1
##
# Also the # of UMIs per cell and the # of Genes per cell are used
mo268[["percent.mt"]] <- PercentageFeatureSet(mo268, pattern = "^MT-")
VlnPlot(mo268,c("nCount_RNA","nFeature_RNA","percent.mt"))

#There is no "explicit" way to select the thresholds

P1<-ggscatterhist(mo268@meta.data, x = "nCount_RNA", y = "nFeature_RNA",color = "percent.mt", 
                  size = 1, alpha = 0.6, 
                  margin.params = list(fill="black",color = "black", size = 0.2))

P2<-ggdensity(mo268@meta.data, x = "percent.mt",
              add = "mean", rug = TRUE)+geom_vline(xintercept = c(5,10,20),linetype="dotted",color = "red")


P1+P2

plot(mo268@meta.data$nCount_RNA,mo268@meta.data$nFeature_RNA,pch=16,bty="n")
abline(h=c(1000,7000),v=c(2500,7e04),col="red",lty=2)

# Filtering based on QC parameters
mo268 <- subset(mo268, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & nCount_RNA > 2500 & nCount_RNA < 5e04 & percent.mt < 10)

# Normalize and Scale Data
mo268<-NormalizeData(mo268)
mo268 <- ScaleData(mo268, features = rownames(mo268))

## Cell Cycle Analysis
mo268<-CellCycleScoring(mo268,s.features = cc.genes$s.genes,g2m.features = cc.genes$g2m.genes)

# Using PCA, we see if there is an effect caussed due to cell cycle phase
mo268<-RunPCA(mo268,features = c(cc.genes$s.genes,cc.genes$g2m.genes))
DimPlot(mo268,group.by = "Phase")
mo268$CC.Difference <- mo268$S.Score - mo268$G2M.Score


# Distribution of FANCA gene

mo268@meta.data$FANCA <- mo268@assays$RNA@data["FANCA",]

mo268@meta.data$Condition<-"Healthy"
mo268@meta.data$Project<-"Healthy"





## Load healthy object

healthy<-readRDS("DATA/seurat_young.rds")

## label transfer and anchoring from healthy donor information

fanconi.anchors <- FindTransferAnchors(reference = healthy, query = mo268, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = fanconi.anchors, refdata = healthy$CellType, 
                            dims = 1:30)
mo268 <- AddMetaData(mo268, metadata = predictions)

mo268$CellType<-mo268$predicted.id

#DimPlot(mo268, group.by = "predicted.id", label=TRUE)
#dev.off()

saveRDS(mo268, "DATA/mo268_names.rds")
