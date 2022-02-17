###############################
#16/10/2020
# DEA in integrated Data
#Miren Lasaga
##############################


#>> Load libraries

library(Seurat)
library(ggplot2)
library(ggpubr)

#>> Work directory and seed

setwd("/datos_2/FANCONI")
set.seed(1234567)

#>> Load integrated data

integration<-readRDS("DATA/integration_new_healthy.RDS")

#>> Re-order the identities of the object

integration@meta.data$CellType<-factor(x=integration@meta.data$CellType, levels=c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils"))

Idents(integration)<-integration$CellType

integration@meta.data$Project[integration@meta.data$orig.ident=="Fanconi"]<-"Fanconi"


# CREATE CORRECTION GROUP
integration$FANCA<-integration@assays$RNA@data["FANCA",]
integration@meta.data$Correction <- as.character(integration@meta.data$FANCA > 0)


integration@meta.data$Correction[integration@meta.data$Correction==TRUE] <- "Corrected"
integration@meta.data$Correction[integration@meta.data$Correction==FALSE] <- "Uncorrected"


integration@meta.data$Correction<- as.factor(integration@meta.data$Correction)

 # Normalize and Scale Data
 integration<-NormalizeData(integration)
 integration <- ScaleData(integration, features = rownames(integration))
 
 ## Cell Cycle Analysis
 integration<-CellCycleScoring(integration,s.features = cc.genes$s.genes,g2m.features = cc.genes$g2m.genes)
 
 # Using PCA, we see if there is an effect caussed due to cell cycle phase
 integration<-RunPCA(integration,features = c(cc.genes$s.genes,cc.genes$g2m.genes))
 DimPlot(integration,group.by = "Phase")
 integration$CC.Difference <- integration$S.Score - integration$G2M.Score
 
 integration <- RunPCA(integration, npcs = 50, verbose = FALSE)
 
  integration<-SCTransform(integration,vars.to.regress = c("S.Score","G2M.Score","percent.mt","nFeature_RNA"))
 
 
 # PCA
 integration <- RunPCA(integration, npcs = 50, verbose = FALSE)
 P1<-DimPlot(fanconi_2008,group.by = "Phase")
 P2<-DimPlot(fanconi_2008)



#>> Run PC analysis

integration <- RunPCA(integration, verbose = FALSE)
integration <- RunUMAP(integration, dims = 1:30)
plots <- DimPlot(integration, group.by = c("Project", "CellType"))

pdf("Integration_updated/Visualizations/UMAP_celltype_project.pdf", width=10)
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 4, byrow = TRUE, 
    override.aes = list(size = 1)))
dev.off()

#>> Check MT %

pdf(file="Mito_percent_UMAP.pdf", width=12, height=6)
p0 = DimPlot(object = integration, reduction = "umap", group.by = "percent.mito")
ggplot(p0$data, aes(p0$data$UMAP_1, p0$data$UMAP_2)) +
    geom_point(aes(colour = integration@meta.data$percent.mt)) +
    scale_colour_gradient2()
dev.off()

#>> Check cell cycle
pdf("INTEGRATION/umap_cellclycle.pdf")
DimPlot(integration,group.by = "Phase")
dev.off()

#>> Clustering in the integrated data


integration@meta.data <- integration@meta.data[, -which(colnames(integration@meta.data) %in% 'seurat_clusters')]

integration@meta.data <- integration@meta.data[, -which(colnames(integration@meta.data) %in% 'Cluster')]

ElbowPlot(integration,ndims = 50)

## UMAP
integration <- RunUMAP(integration, reduction = "pca", dims = 1:22)

integration <- FindNeighbors(integration, reduction = "pca", dims = 1:23)
integration <- FindClusters(integration, resolution = 0.8)

#>> Transform Numbers to Letters
integration$Cluster<-LETTERS[as.numeric(as.vector(integration$seurat_clusters))+1]
LX<-length(levels(integration$seurat_clusters))
integration$Cluster<-factor(integration$Cluster,levels = LETTERS[1:LX])

DefaultAssay(integration)<-"RNA"


library(ComplexHeatmap)
library(RColorBrewer)

my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")
pdf("umap_celltype.pdf")
DimPlot(integration, group.by="CellType", label=TRUE, cols=c(my_palette1[2],my_palette2[1],
                       my_palette2[2],my_palette2[3],my_palette2[4],
                       my_palette2[6],my_palette2[7],
                       my_palette2[8], my_palette2[9],my_palette2[10],my_palette1[4],my_palette2[5]), size=7)
dev.off()


pdf("umap_condition.pdf")
DimPlot(integration, order=c("Corrected","Uncorrected","Healthy"),group.by="Condition", pt.size=1, cols=c("#FFD700","#00BFC4","#F8766D"))#FFD700"
dev.off()


# plot fanca in a boxplot only of the corrected cells

fanca_filtered<-subset(integration, subset=FANCA>0)



pdf("FANCA_expression_filtered.pdf", width=10)
p<-ggplot(fanca_filtered@meta.data, aes(x=CellType, y=FANCA, fill=Sample)) + 
    geom_boxplot()+
    scale_fill_manual(values=c("#F8766D","#FFD700"))+
    labs(title="Expression of FANCA by cell type",x = "Cell_type", y = "FANCA")

p +theme_bw()
dev.off()
geom_jitter(shape=16, position=position_jitter(0.2))
#> Plot number of cells and percentage by cell type

## plot healthy expressing and not expressing FANCA

sano<-subset(integration, subset=Project=="Healthy")

sano$FANCA_expression<-sano$FANCA>0
sano@meta.data$FANCA_expression[sano@meta.data$FANCA_expression==TRUE] <- "Fanca_expressed"
sano@meta.data$FANCA_expression[sano@meta.data$FANCA_expression==FALSE] <- "Fanca_not_expressed"


sano@meta.data$FANCA_expression<- as.factor(sano@meta.data$FANCA_expression)


for (i in 1:dim(integration)[1]){
if (is.na(integration$Condition[i]))
{integration$Condition[i]<-as.factor(integration$Correction[i])}
}

pdf("Fanconi_2008/integration/fanca_expression_by_project.pdf",width=8)

p1<- ggplot(integration@meta.data, aes(x=CellType, fill=Project)) + 
  geom_bar(position="fill")+
  scale_fill_manual(values=c("#9400D3","#FFD700"))+#9400D3
    geom_text(aes(label=..count..),stat='count',position=position_fill(vjust=0.5),size=3)+
  labs(title="% of correction per Cell Type",x = "Cell_type", y = "%")
p1+coord_flip()+theme(text = element_text(size=12))
dev.off()



## plot celltype by sample
library(RColorBrewer)
colors <- c(brewer.pal(n=8, name = 'Dark2'), brewer.pal(n=8, name='Set3'),brewer.pal(n=12, name = 'Paired'))

color_types <- colors[1:length(levels(as.factor(integration$CellType)))]
names(color_types) <- levels(integration$CellType)

p1<-ggplot(integration@meta.data, aes(x=Combination, fill=CellType)) + geom_bar(position="fill") +
scale_fill_manual(values=c(my_palette1[2],my_palette2[1],
                       my_palette2[2],my_palette2[3],my_palette2[4],
                       my_palette2[6],my_palette2[7],
                       my_palette2[8], my_palette2[9],my_palette2[10],my_palette1[4],my_palette2[5]))+
  labs(title="% of cells types by Individual",x = "Individual", y = "%")
 
ggsave(plot = p1, filename = "Integration_updated/Visualizations/percentage_by_cell_type.pdf", width=10, height=10, dpi=30)


#>> Find markers healthy vs fanconi patient by cellType
setwd("/datos_2/FANCONI/Integration_new_healthy")
library(biomaRt)
library(ggplot2)
library(clusterProfiler)
library('org.Hs.eg.db') 
set.seed(1234567)

Cell_type<-c("HSC","LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils") #HSC cab't be included, any corrected cell

for(i in 1:length(Cell_type)) {
  
disease_markers<-FindMarkers(integration, ident.1="Healthy", ident.2="Uncorrected", group.by= "Condition", logfc.threshold = 0, subset.ident = Cell_type[i])

sig_genes_up <- which(disease_markers$p_val_adj<0.05 & disease_markers$avg_logFC>0.25)
sig_genes_down <- which(disease_markers$p_val_adj<0.05 & disease_markers$avg_logFC< -0.25)
disease_markers$sig_genes <- rep(0,nrow(disease_markers))
  disease_markers$sig_genes[disease_markers$avg_logFC>0] <-0.5
  disease_markers$sig_genes[disease_markers$avg_logFC<0] <- -0.5
  disease_markers$sig_genes[sig_genes_up] <- 1
  disease_markers$sig_genes[sig_genes_down] <- -1
  colnames(disease_markers)<-paste(colnames(disease_markers), Cell_type[i], sep = "_")
write.table(disease_markers, paste0("Healthy_Uncorrected_",Cell_type[i],".txt", sep=""),sep="\t",quote=F,row.names=T)

}
## MERGE ALL THE DATA FROM DEA AND PLOT HEATMAP

## MERGE ALL THE DATA FROM DEA AND PLOT HEATMAP




dea_HSC_fanconi<-read.table("Healthy_Uncorrected_HSC.txt", sep="\t", header=T)
dea_LMPP_fanconi<-read.table("Healthy_Uncorrected_LMPP.txt", sep="\t", header=T)
dea_Cycling_LMPP_fanconi<-read.table("Healthy_Uncorrected_Cycling_LMPP.txt", sep="\t", header=T)
dea_GMP1_fanconi<-read.table("Healthy_Uncorrected_GMP1.txt", sep="\t", header=T)
dea_GMP2_fanconi<-read.table("Healthy_Uncorrected_GMP2.txt", sep="\t", header=T)
dea_Monocytes_fanconi<-read.table("Healthy_Uncorrected_Monocytes.txt", sep="\t", header=T)
dea_DC_fanconi<-read.table("Healthy_Uncorrected_DC.txt", sep="\t", header=T)
dea_CLP_fanconi<-read.table("Healthy_Uncorrected_CLP.txt", sep="\t", header=T)
dea_PreB_fanconi<-read.table("Healthy_Uncorrected_PreB.txt", sep="\t", header=T)
dea_MEP_fanconi<-read.table("Healthy_Uncorrected_MEP.txt", sep="\t", header=T)
dea_Erythroid_fanconi<-read.table("Healthy_Uncorrected_Erythroid.txt", sep="\t", header=T)
dea_Basophils_fanconi<-read.table("Healthy_Uncorrected_Basophils.txt", sep="\t", header=T)

 

MyMerge       <- function(x, y){
  df            <- merge(x, y, by= "row.names", all.x= T, all.y= T)
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)
}
dat           <- Reduce(MyMerge, list(dea_HSC_fanconi,dea_LMPP_fanconi,dea_Cycling_LMPP_fanconi, dea_GMP1_fanconi, dea_GMP2_fanconi,dea_Monocytes_fanconi,dea_DC_fanconi,dea_CLP_fanconi,dea_PreB_fanconi,dea_MEP_fanconi,dea_Erythroid_fanconi,dea_Basophils_fanconi))

sig_genes<-read.table("union_3_samples_sig_genes.txt",sep="\t")
sig<-as.character(sig_genes$x)

fanconi_2006 <- dat[,grep("sig_genes",colnames(dat))]
colnames(fanconi_2006) <- gsub("sig_genes_","",colnames(fanconi_2006))
#data_old<-data_old[intersection,]
 fanconi_2006[is.na(fanconi_2006)] <- 0


uncorrected_sig<-uncorrected_2006[sig,]
corrected_sig<-corrected_2006[sig,]
fanconi_sig<-fanconi_2006[sig,]

data_2006<-cbind(fanconi_sig, uncorrected_sig, corrected_sig)

write.table(data_2004, "Fanconi_eltrombopag/merge_3_contrast.txt", sep="\t", row.names=TRUE)


merge<-cbind(data_2008,data_2006,data_2004)

merge[is.na(merge)] <- 0

metadata_all <- data.frame(ID=colnames(all),
                         cell_type=factor(colnames(all),levels=c("HSC","LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils")),
                          group=factor(c(rep(c("corrected_vs_uncorrected","healthy_vs_uncorrected","healthy_vs_corrected"), each=10),rep(c("corrected_vs_uncorrected","healthy_vs_uncorrected","healthy_vs_corrected"), each=11),rep(c("corrected_vs_uncorrected","healthy_vs_uncorrected","healthy_vs_corrected"), each=11),rep(c("corrected_vs_uncorrected","healthy_vs_uncorrected","healthy_vs_corrected"),each=12)), levels=c("corrected_vs_uncorrected","healthy_vs_uncorrected","healthy_vs_corrected")),
                          sample=factor(c(rep(c("sample_2008"),each=30),rep(c("sample_2006"),each=33),rep(c("sample_2004"),each=33),rep(c("sample_2002"),each=36)),level=c("sample_2002","sample_2004","sample_2006","sample_2008")))


ha1 <- HeatmapAnnotation(
  Cell_type = metadata_all$cell_type,
  Group=metadata_all$group,
  Sample=metadata_all$sample,
  col = list(Cell_type = c("HSC"=my_palette1[2],"LMPP"=my_palette2[1],
                       "Cycling_LMPP"=my_palette2[2],"GMP1"=my_palette2[3],"GMP2"=my_palette2[4],
                       "Monocytes"=my_palette2[6],"DC"=my_palette2[7],
                       "CLP"=my_palette2[8], "PreB"=my_palette2[9],"MEP"=my_palette2[10],"Erythroid"=my_palette1[4],"Basophils"=my_palette2[5]),
              Group= c("corrected_vs_uncorrected"="#CD853F" ,"healthy_vs_uncorrected"= "#8B008B", "healthy_vs_corrected"="#6B8E23" ),
              Sample= c("sample_2008"="#FF6347", "sample_2006"="#6495ED", "sample_2004"="#008080", "sample_2002"="lightpink")),
  show_annotation_name = TRUE
)
#plot
#jpeg("DEA/MM_C/heatmap_binary_sig_genes_MM_Control.jpg", width=5000, height=5000, res=600)

## show target gene names in heatmap

target_TF <- c("DNMT1","MCM2","PCNA","MCM10","MCM7","UBE2T","DKC1","POLR3K","PARP1","FANCI","RAD51C","MSH6")
target_TF_position <- vector()
for(i in 1:length(target_TF)){
  target_TF_position <- c(target_TF_position,grep(paste("_",target_TF[i],sep=""),rownames(merge)))
}
rha = rowAnnotation(foo = anno_mark(at = target_TF_position,
                                    labels = rownames(merge)[target_TF_position],
                                    labels_gp = gpar(fontsize = 6)))


pdf("INTERSECTION_sig_3_samples.pdf")
#cell_type > group > sample
Heatmap(all, name = "DEGs genes",
        col = structure(c(brewer.pal(n = 8, name = "YlGnBu")[c(7,5)], "gray", brewer.pal(n = 8, name = "YlOrRd")[c(4,7)]), names=c("-1","-0.5","0","0.5","1")),
        #col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 3), 
        cluster_rows = TRUE, cluster_columns = FALSE, column_order=order(metadata_all$group,metadata_all$cell_type, metadata_all$sample) , use_raster=TRUE) 
dev.off()

#cell_type > group > sample
Heatmap(def, name = "DEGs genes",
        col = structure(c(brewer.pal(n = 8, name = "YlGnBu")[c(7,5)], "gray", brewer.pal(n = 8, name = "YlOrRd")[c(4,7)]), names=c("-1","-0.5","0","0.5","1")),
        #col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = TRUE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 3), 
        cluster_rows = TRUE, cluster_columns = FALSE, column_order=order(metadata_all$group, metadata_all$cell_type, metadata_all$sample) , use_raster=TRUE) 


## take all the genes that are significant in at least 1 celltype

background <- unique(c(dea_LMPP_fanconi, dea_Cycling_LMPP_fanconi, dea_GMP1_fanconi, dea_GMP2_fanconi,dea_Monocytes_fanconi,dea_DC_fanconi,dea_CLP_fanconi,dea_PreB_fanconi,dea_MEP_fanconi,dea_Erythroid_fanconi,dea_Basophils_fanconi))

sig_genes<-integration@assays$RNA@data[background,]

design_matrix<-integration[background,]

library(ComplexHeatmap)
library(RColorBrewer)

my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Cell_type = metadata$CellType,
  signature = metadata$signature,
  col = list(signature = c("FANCA_expressing_cell"=my_palette1[3],"Signature_positive_cell"=my_palette1[1],"Signature_negative_cell"=my_palette1[7]),
             Cell_type = c("HSC"=my_palette1[2],"LMPP"=my_palette2[1],
                       "Cycling_LMPP"=my_palette2[2],"GMP1"=my_palette2[3],"GMP2"=my_palette2[4],
                       "Monocytes"=my_palette2[6],"DC"=my_palette2[7],
                       "CLP"=my_palette2[8], "PreB"=my_palette2[9],"MEP"=my_palette2[10],"Erythroid"=my_palette1[4],"Basophils"=my_palette2[5])),
  show_annotation_name = TRUE
  

)

sig_genes2<-as.matrix(sig_genes)
data_to_plot2<- t(scale(t(data_to_plot2)))

#plot
pdf("sig_genes_expression_corrected.pdf") 
#jpeg("Results/MM/expression_heatmap_genotype.jpg", width=5000, height=5000, res=600)#

Heatmap(def, name = "DEGs corrected vs healthy",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 1), 
        cluster_rows = TRUE, cluster_columns = FALSE,column_order=order(metadata$signature, metadata$CellType), use_raster=TRUE)



dev.off()


