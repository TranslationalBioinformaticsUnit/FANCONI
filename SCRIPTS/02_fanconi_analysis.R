###### analysis of the fanconi sample

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
fanconi<-Read10X("DATA/fanconi_eltrombopag")

 fanconi <- CreateSeuratObject(counts = fanconi, project = "Fanconi", min.cells = 3, min.features = 400)
 
 # Initial QC Netrics
 # We get the % of UMIs mapped to mitochondrial genes as a way to identify "broken/stressed/ cells
 ## Reference
 #Ilicic, T., Kim, J.K., Kolodziejczyk, A.A. et al. Classification of low quality cells from single-cell RNA-seq data. Genome Biol 17, 29 (2016). 
 #https://doi.org/10.1186/s13059-016-0888-1
 ##
 # Also the # of UMIs per cell and the # of Genes per cell are used
 fanconi[["percent.mt"]] <- PercentageFeatureSet(fanconi, pattern = "^MT-")
 VlnPlot(fanconi,c("nCount_RNA","nFeature_RNA","percent.mt"))
 
 #There is no "explicit" way to select the thresholds
 
 P1<-ggscatterhist(fanconi@meta.data, x = "nCount_RNA", y = "nFeature_RNA",color = "percent.mt", 
                   size = 1, alpha = 0.6, 
                   margin.params = list(fill="black",color = "black", size = 0.2))
 
 P2<-ggdensity(fanconi@meta.data, x = "percent.mt",
               add = "mean", rug = TRUE)+geom_vline(xintercept = c(5,10,20),linetype="dotted",color = "red")
 
 
 P1+P2
 
 plot(fanconi@meta.data$nCount_RNA,fanconi@meta.data$nFeature_RNA,pch=16,bty="n")
 abline(h=c(1000,7000),v=c(2500,7e04),col="red",lty=2)
 
 # Filtering based on QC parameters
 fanconi <- subset(fanconi, subset = nFeature_RNA > 700 & nFeature_RNA < 7000 & nCount_RNA > 1000 & nCount_RNA < 5e04 & percent.mt < 10)  # usually we filter nFeatures and nCouts 10 % up and down and MT at 10%
 
 # Normalize and Scale Data
 fanconi<-NormalizeData(fanconi)
 fanconi <- ScaleData(fanconi, features = rownames(fanconi))
 
 ## Cell Cycle Analysis
 fanconi<-CellCycleScoring(fanconi,s.features = cc.genes$s.genes,g2m.features = cc.genes$g2m.genes)
 
 # Using PCA, we see if there is an effect caussed due to cell cycle phase
 fanconi<-RunPCA(fanconi,features = c(cc.genes$s.genes,cc.genes$g2m.genes))
 DimPlot(fanconi,group.by = "Phase")
 fanconi$CC.Difference <- fanconi$S.Score - fanconi$G2M.Score

 
 # SC Transform (https://satijalab.org/seurat/v3.1/sctransform_vignette.html)
 fanconi<-SCTransform(fanconi,vars.to.regress = c("S.Score","G2M.Score","percent.mt","nFeature_RNA"))
 
 
 # PCA
 fanconi <- RunPCA(fanconi, npcs = 50, verbose = FALSE)
 P1<-DimPlot(fanconi,group.by = "Phase")
 P2<-DimPlot(fanconi)
 P1+P2
 
 # Select components to use in further steps
 ElbowPlot(fanconi,ndims = 50)
 
 ## UMAP
 fanconi <- RunUMAP(fanconi, reduction = "pca", dims = 1:16)
 
 ## Clustering
 fanconi <- FindNeighbors(fanconi, reduction = "pca", dims = 1:12)
 fanconi <- FindClusters(fanconi, resolution = 0.8, n.start = 1000) ## n.start=1000 para obtener resultados similares dependiendo del orden de células
 
 ## Label transfer
 healthy<-readRDS("DATA/seurat_young.rds")

## label transfer and anchoring from healthy donor information

fanconi.anchors <- FindTransferAnchors(reference = healthy, query = fanconi, 
                                        dims = 1:30,project.query = T)
predictions <- TransferData(anchorset = fanconi.anchors, refdata = healthy$CellType, 
                            dims = 1:30)
fanconi <- AddMetaData(fanconi, metadata = predictions)

# MT %
pdf("RESULTS/clustering_Cell_type_prediction_MT.pdf")
p0 = DimPlot(object = fanconi, reduction = "umap", group.by = "percent.mt")
ggplot(p0$data, aes(p0$data$UMAP_1, p0$data$UMAP_2)) +
  geom_point(aes(colour = fanconi@meta.data$percent.mt)) +
  scale_colour_gradient2()
dev.off()


#fanconi<-readRDS("/datos_2/FANCONI/DATA/Fanconi_CD34_MO.rds")

fanconi@meta.data$CellType<- fanconi@meta.data$predicted.id


fanconi@meta.data$CellType<-factor(x=fanconi@meta.data$CellType, levels=c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils"))

fanconi<-SetIdent(fanconi, value="CellType")

levels(x = fanconi)<-c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils")






# CREATE CORRECTION GROUP
fanconi$FANCA<-fanconi@assays$RNA@data["FANCA",]
fanconi@meta.data$Correction <- as.character(fanconi@meta.data$FANCA > 0)


fanconi@meta.data$Correction[fanconi@meta.data$Correction==TRUE] <- "Corrected"
fanconi@meta.data$Correction[fanconi@meta.data$Correction==FALSE] <- "Uncorrected"


fanconi@meta.data$Correction<- as.factor(fanconi@meta.data$Correction)


# PLOT UMAP BY CELLTYPE AND CORRECTION

pdf("Fanconi_eltrombopag/umap_celltype.pdf")
DimPlot(fanconi, group.by="CellType", label=TRUE, cols=c(my_palette1[2],my_palette2[1],
                       my_palette2[2],my_palette2[3],my_palette2[4],
                       my_palette2[6],my_palette2[7],
                       my_palette2[8], my_palette2[9],my_palette2[10],my_palette1[4],my_palette2[5]), size=10)
dev.off()



pdf("Fanconi_eltrombopag/umap_correction.pdf")
DimPlot(fanconi, group.by="Correction")
dev.off()

# barplots with % of correction by celltype


p1<- ggplot(fanconi@meta.data, aes(x=CellType, fill=Correction)) + 
  geom_bar(position="fill")+
  #scale_fill_manual(values=c("#D95F02","#7570B3"))+
    geom_text(aes(label=..count..),stat='count',position=position_fill(vjust=0.5),size=7)+
  labs(title="% of correction per Cell Type",x = "Cell_type", y = "%")
pdf("Fanconi_eltrombopag/percentage_cell_type.pdf")
p1+coord_flip()+theme(axis.text=element_text(size=10))
dev.off()



p1<-ggplot(fanconi@meta.data, aes(x=CellType, fill=Correction)) +
geom_bar() + 
labs(title="Number of cells per Cell Type",x = "Cell_type", y = "Number of cells")+
stat_count(aes(label = ..count..), geom = "text")
ggsave(plot = p1, height = 7, width = 9, dpi=600, filename = "PLOTS/cells_by_cell_type.pdf")



# plot fanca in a boxplot only of the corrected cells

corrected<-subset(fanconi, subset=Correction=="Corrected")

p<-ggplot(corrected@meta.data, aes(x=CellType, y=FANCA, fill=Sample)) + 
    geom_boxplot()+
    #scale_fill_manual(values=c("#D95F02"))+
    labs(title="Expression of FANCA by cell type",x = "Cell_type", y = "FANCA")

pdf("Fanconi_eltrombopag/FANCA_expression.pdf", width=10)
p + geom_jitter(shape=16, position=position_jitter(0.2))+theme_bw()
dev.off()


#Differential expression analysis with all the genes logFC>0 per cell type  


Cell_type<-c("HSC","LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","MEP","Erythroid","Basophils") #HSC cab't be included, any corrected cell

for(i in 1:length(Cell_type)) {
  
disease_markers<-FindMarkers(fanconi, ident.1="Corrected", ident.2="Uncorrected", group.by= "Correction", logfc.threshold = 0, subset.ident = Cell_type[i])
sig_genes_up <- which(disease_markers$p_val_adj<0.05 & disease_markers$avg_logFC>0.25)
sig_genes_down <- which(disease_markers$p_val_adj<0.05 & disease_markers$avg_logFC< -0.25)
disease_markers$sig_genes <- rep(0,nrow(disease_markers))
  disease_markers$sig_genes[disease_markers$avg_logFC>0] <-0.5
  disease_markers$sig_genes[disease_markers$avg_logFC<0] <- -0.5
  disease_markers$sig_genes[sig_genes_up] <- 1
  disease_markers$sig_genes[sig_genes_down] <- -1
  colnames(disease_markers)<-paste(colnames(disease_markers), Cell_type[i], sep = "_")
write.table(disease_markers, paste0("Fanconi_eltrombopag/Corrected_vs_uncorrected_",Cell_type[i],".txt", sep=""),sep="\t",quote=F,row.names=T)
markers<-NULL
}

## MERGE ALL THE DATA FROM DEA AND PLOT HEATMAP
#dea_HSC_fanconi<-read.table("Fanconi_eltrombopag/Corrected_vs_uncorrected_HSC.txt", sep="\t", header=T)
dea_LMPP_fanconi<-read.table("Fanconi_sample/DEA_all/Corrected_vs_uncorrected_LMPP.txt", sep="\t", header=T)
dea_Cycling_LMPP_fanconi<-read.table("Fanconi_sample/DEA_all/Corrected_vs_uncorrected_Cycling_LMPP.txt", sep="\t", header=T)
dea_GMP1_fanconi<-read.table("Fanconi_sample/DEA_all/Corrected_vs_uncorrected_GMP1.txt", sep="\t", header=T)
dea_GMP2_fanconi<-read.table("Fanconi_sample/DEA_all/Corrected_vs_uncorrected_GMP2.txt", sep="\t", header=T)
dea_Monocytes_fanconi<-read.table("Fanconi_sample/DEA_all/Corrected_vs_uncorrected_Monocytes.txt", sep="\t", header=T)
dea_DC_fanconi<-read.table("Fanconi_sample/DEA_all/Corrected_vs_uncorrected_DC.txt", sep="\t", header=T)
dea_CLP_fanconi<-read.table("Fanconi_sample/DEA_all/Corrected_vs_uncorrected_CLP.txt", sep="\t", header=T)
dea_PreB_fanconi<-read.table("Fanconi_sample/DEA_all/Corrected_vs_uncorrected_PreB.txt", sep="\t", header=T)
dea_MEP_fanconi<-read.table("Fanconi_sample/DEA_all/Corrected_vs_uncorrected_MEP.txt", sep="\t", header=T)
dea_Erythroid_fanconi<-read.table("Fanconi_sample/DEA_all/Corrected_vs_uncorrected_Erythroid.txt", sep="\t", header=T)
dea_Basophils_fanconi<-read.table("Fanconi_sample/DEA_all/Corrected_vs_uncorrected_Basophils.txt", sep="\t", header=T)


MyMerge       <- function(x, y){
  df            <- merge(x, y, by= "row.names", all.x= T, all.y= T)
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)
}
dat           <- Reduce(MyMerge, list( dea_LMPP_fanconi, dea_Cycling_LMPP_fanconi, dea_GMP1_fanconi, dea_GMP2_fanconi,dea_Monocytes_fanconi,dea_DC_fanconi,dea_CLP_fanconi, dea_PreB_fanconi,dea_MEP_fanconi,dea_Erythroid_fanconi,dea_Basophils_fanconi))


## in dat we merge all the ifnormation of all the cell types

# union or intersection is the list of gene of interest!!

data_fanconi<-dat[union,]

data_to_plot_fanconi<- data_fanconi[,grep("sig_genes",colnames(data_fanconi))]
colnames(data_to_plot_fanconi) <- gsub("sig_genes_","",colnames(data_to_plot_fanconi))
#data_old<-data_old[intersection,]
 data_to_plot_fanconi[is.na(data_to_plot_fanconi)] <- 0

ha1 <- HeatmapAnnotation(
  Cell_type = colnames(data_to_plot_fanconi),
  col = list(Cell_type = c("HSC"=my_palette1[2],"LMPP"=my_palette2[1],
                       "Cycling_LMPP"=my_palette2[2],"GMP1"=my_palette2[3],"GMP2"=my_palette2[4],
                       "Monocytes"=my_palette2[6],"DC"=my_palette2[7],
                       "CLP"=my_palette2[8], "PreB"=my_palette2[9],"MEP"=my_palette2[10],"Erythroid"=my_palette1[4],"Basophils"=my_palette2[5])),
  show_annotation_name = TRUE
)
library(dendextend)
row_dend = as.dendrogram(hclust(dist(data_to_plot_fanconi)))
col_dend = as.dendrogram(hclust(dist(t(data_to_plot_fanconi))))
#plot
#jpeg("DEA/MM_C/heatmap_binary_sig_genes_MM_Control.jpg", width=5000, height=5000, res=600)
pdf("/datos_2/FANCONI/Fanconi_sample/PLOTS/heatmap_binary_sig_genes_ordered.pdf", width = 8, height = 8)
Heatmap(data_to_plot_fanconi, name = "DEGs between corrected and uncorrected",
        col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, use_raster=TRUE)
dev.off()


