#----- Plot boxplot of each sample in public data. target genes: MYC pathway genes

setwd("/datos_2/FANCONI/Public_data/Analysis/MYC_plots")
library(Seurat)
library(ggplot2)

#------Define target genes

target_genes<-read.table("/datos_2/FANCONI/AUCell_analysis/MYC_pathway/genes_pathway.txt")
target_genes<-as.character(target_genes$x)

genes_to_remove<-c("PPP2R4")

target_genes<-target_genes[!(target_genes %in% genes_to_remove)]

#------ Load public data seurat object

data<-readRDS("/datos_2/FANCONI/Public_data/integrated_cellnames.RDS")


#par(mfrow=c(1,1))

num_columns=1
num_rows=1
num_pages=1
val=num_columns*num_rows*num_pages

vplots=list()
for(i in 1:length(target_genes)){

  integration$target_genes<-integration@assays$RNA@data[target_genes[i],]
  
   #data2<-subset(integration, subset=target_genes>0)
   
    p<- ggplot(data = integration@meta.data, aes(x=Sample, y=target_genes, fill=Sample)) + 
         geom_boxplot() +
    labs(x = "Cell_type", y = target_genes[i])+
    ggtitle(paste(target_genes[i]," - normalized gene expression", sep=""))
    
   
    vplots[[i]]=p 
    
 }



pdf("Analysis/Visualizations/MYC_pathway_genes_boxplot.pdf")
marrangeGrob(vplots, nrow = 1, ncol = 1)
dev.off()

#plot MYC expression
data$MYC<-data@assays$RNA@data["MYC",]
  
   data2<-subset(data, subset=MYC>0)

 pdf("MYC_expression_sample.pdf")  
    p<- ggplot(data = data@meta.data, aes(x=Sample, y=MYC, fill=Sample)) + 
         geom_boxplot() +
    labs(x = "Cell_type", y = "MYC")+
    ggtitle("MYC - normalized gene expression")
 dev.off()  



### load our contrast for monocytes and erythroids our data and public data

dea_Monocytes_2006<-read.table("Fanconi_sample/DEA_all/Corrected_vs_uncorrected_Monocytes.txt", sep="\t", header=T)
dea_Erythroids_2006<-read.table("Fanconi_sample/DEA_all/Corrected_vs_uncorrected_Erythroid.txt", sep="\t", header=T)
dea_Monocytes_2004<-read.table("Fanconi_eltrombopag/Corrected_vs_uncorrected_Monocytes.txt", sep="\t", header=T)
dea_Erythroids_2004<-read.table("Fanconi_eltrombopag/Corrected_vs_uncorrected_Erythroid.txt", sep="\t", header=T)
dea_Monocytes_2008<-read.table("Fanconi_2008/Corrected_vs_uncorrected_Monocytes.txt", sep="\t", header=T)
dea_Erythroids_2008<-read.table("Fanconi_2008/Corrected_vs_uncorrected_Erythroid.txt", sep="\t", header=T)
dea_Monocytes_public<-read.table("Public_data/Analysis/DEA/Healthy_fanconi_Monocytes.txt", sep="\t", header=T)
dea_Erythroids_public<-read.table("Public_data/Analysis/DEA/Healthy_fanconi_Erythroid.txt", sep="\t", header=T)

sig_genes<-read.table("union_3_samples_sig_genes.txt",sep="\t")
sig<-as.character(sig_genes$x)



# merge data of each sample and take info of sig genes

fanconi_2006<-merge(dea_Monocytes_2006, dea_Erythroids_2006, by.x=0, by.y=0, all.x=TRUE, all.y=TRUE)
rownames(fanconi_2006)<-fanconi_2006$Row.names
fanconi_2006[is.na(fanconi_2006)] <- 0
fanconi_2006 <- fanconi_2006[,grep("sig_genes",colnames(fanconi_2006))]
colnames(fanconi_2006) <- gsub("sig_genes_","",colnames(fanconi_2006))
fanconi_2006_sig<-fanconi_2006[union,]
 
fanconi_2004<-merge(dea_Monocytes_2004, dea_Erythroids_2004, by.x=0, by.y=0, all.x=TRUE, all.y=TRUE)
rownames(fanconi_2004)<-fanconi_2004$Row.names
fanconi_2004[is.na(fanconi_2004)] <- 0
fanconi_2004 <- fanconi_2004[,grep("sig_genes",colnames(fanconi_2004))]
colnames(fanconi_2004) <- gsub("sig_genes_","",colnames(fanconi_2004))
fanconi_2004_sig<-fanconi_2004[union,]


fanconi_2008<-merge(dea_Monocytes_2008, dea_Erythroids_2008, by.x=0, by.y=0, all.x=TRUE, all.y=TRUE)
rownames(fanconi_2008)<-fanconi_2008$Row.names
fanconi_2008[is.na(fanconi_2008)] <- 0
fanconi_2008 <- fanconi_2008[,grep("sig_genes",colnames(fanconi_2008))]
colnames(fanconi_2008) <- gsub("sig_genes_","",colnames(fanconi_2008))
fanconi_2008_sig<-fanconi_2008[union,]


public_data<-merge(dea_Monocytes_public, dea_Erythroids_public, by.x=0, by.y=0, all.x=TRUE, all.y=TRUE)
rownames(public_data)<-public_data$Row.names
public_data[is.na(public_data)] <- 0
public_data <- public_data[,grep("sig_genes",colnames(public_data))]
colnames(public_data) <- gsub("sig_genes_","",colnames(public_data))
public_data_sig<-public_data[union,]





merge<-cbind(fanconi_2004_sig,fanconi_2006_sig,fanconi_2008_sig, public_data_sig)

## merge 2 with the intersection

sig_genes<-read.table("intersection_sig_genes_3_samples.txt",sep="\t")
sig<-as.character(sig_genes$x)

merge2<-merge[sig,]

metadata_all <- data.frame(ID=colnames(merge),
                         cell_type=factor(colnames(merge),levels=c("Monocytes","Erythroid")),
                         sample=factor(c(rep(c("sample_2008"),each=2),rep(c("sample_2006"),each=2),rep(c("sample_2004"),each=2), rep(c("public_data"),each=2)),levels=c("sample_2004","sample_2006","sample_2008", "public_data")))
                         
library(ComplexHeatmap)
library(RColorBrewer)



my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")



ha1 <- HeatmapAnnotation(
  Cell_type = metadata_all$cell_type,
  Sample=metadata_all$sample,
  col = list(Cell_type = c("Monocytes"=my_palette2[6],"Erythroid"=my_palette1[1]),
              Sample= c("sample_2008"="#FF6347", "sample_2006"="#6495ED", "sample_2004"="#008080", "public_data"="#CD853F")),
  show_annotation_name = TRUE
)
#plot
#jpeg("DEA/MM_C/heatmap_binary_sig_genes_MM_Control.jpg", width=5000, height=5000, res=600)

pdf("Public_data/Analysis/DEA/intersection_fanconi_heatmap_order_sample.pdf")
#cell_type > group > sample
Heatmap(merge2, name = "DEGs genes",
        col = structure(c(brewer.pal(n = 8, name = "YlGnBu")[c(7,5)], "gray", brewer.pal(n = 8, name = "YlOrRd")[c(4,7)]), names=c("-1","-0.5","0","0.5","1")),
        #col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = TRUE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 3), 
        cluster_rows = TRUE, cluster_columns = FALSE, column_order=order(metadata_all$sample,metadata_all$cell_type) , use_raster=TRUE) 
dev.off()


## take the DEA results and get the signature
dea_HSC_fanconi<-read.table("Analysis/DEA/Healthy_fanconi_HSC.txt", sep="\t", header=T)
dea_LMPP_fanconi<-read.table("Analysis/DEA/Healthy_fanconi_LMPP.txt", sep="\t", header=T)
dea_Cycling_LMPP_fanconi<-read.table("Analysis/DEA/Healthy_fanconi_Cycling_LMPP.txt", sep="\t", header=T)
dea_GMP1_fanconi<-read.table("Analysis/DEA/Healthy_fanconi_GMP1.txt", sep="\t", header=T)
dea_GMP2_fanconi<-read.table("Analysis/DEA/Healthy_fanconi_GMP2.txt", sep="\t", header=T)
dea_Monocytes_fanconi<-read.table("Analysis/DEA/Healthy_fanconi_Monocytes.txt", sep="\t", header=T)
dea_DC_fanconi<-read.table("Analysis/DEA/Healthy_fanconi_DC.txt", sep="\t", header=T)
dea_CLP_fanconi<-read.table("Analysis/DEA/Healthy_fanconi_CLP.txt", sep="\t", header=T)
#dea_PreB_fanconi<-read.table("Analysis/DEA/Healthy_fanconi_PreB.txt", sep="\t", header=T)
dea_MEP_fanconi<-read.table("Analysis/DEA/Healthy_fanconi_MEP.txt", sep="\t", header=T)
dea_Erythroid_fanconi<-read.table("Analysis/DEA/Healthy_fanconi_Erythroid.txt", sep="\t", header=T)
#dea_Basophils_fanconi<-read.table("Analysis/DEA/Healthy_fanconi_Basophils.txt", sep="\t", header=T)

MyMerge       <- function(x, y){
  df            <- merge(x, y, by= "row.names", all.x= T, all.y= T)
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)
}
dat           <- Reduce(MyMerge, list( dea_HSC_fanconi, dea_LMPP_fanconi, dea_Cycling_LMPP_fanconi, dea_GMP1_fanconi, dea_GMP2_fanconi,dea_Monocytes_fanconi,dea_DC_fanconi,dea_CLP_fanconi,dea_MEP_fanconi,dea_Erythroid_fanconi))

public_data<- dat[,grep("sig_genes",colnames(dat))]
colnames(public_data) <- gsub("sig_genes_","",colnames(public_data))

 public_data[is.na(public_data)] <- 0
 
write.table(public_data, "Analysis/DEA/all_DEA_significance_data.txt", sep="\t")
#

## take only the significant genes in each case
HSC<-rownames(dea_HSC_fanconi[abs(dea_HSC_fanconi$sig_genes)==1,])
LMPP<-rownames(dea_LMPP_fanconi[abs(dea_LMPP_fanconi$sig_genes)==1,])
Cycling_LMPP<-rownames(dea_Cycling_LMPP_fanconi[abs(dea_Cycling_LMPP_fanconi$sig_genes)==1,])
GMP1<-rownames(dea_GMP1_fanconi[abs(dea_GMP1_fanconi$sig_genes)==1,])
GMP2<-rownames(dea_GMP2_fanconi[abs(dea_GMP2_fanconi$sig_genes)==1,])
Monocytes<-rownames(dea_Monocytes_fanconi[abs(dea_Monocytes_fanconi$sig_genes)==1,])
DC<-rownames(dea_DC_fanconi[abs(dea_DC_fanconi$sig_genes)==1,])
CLP<-rownames(dea_CLP_fanconi[abs(dea_CLP_fanconi$sig_genes)==1,])

MEP<-rownames(dea_MEP_fanconi[abs(dea_MEP_fanconi$sig_genes)==1,])
Erythroid<-rownames(dea_Erythroid_fanconi[abs(dea_Erythroid_fanconi$sig_genes)==1,])
#Basophils<-rownames(dea_Basophils_fanconi[abs(dea_Basophils_fanconi$sig_genes)==1,])

sig_genes<-unique(c(HSC,LMPP,Cycling_LMPP,GMP1,GMP2,Monocytes,DC,CLP,MEP,Erythroid))

write.table(sig_genes, "Analysis/DEA/union_sig_genes.txt", sep="\t")



## PLOT heatmap with the signigicant genes

public_data_2<-public_data[sig_genes,]

ha1 <- HeatmapAnnotation(
  Cell_type = colnames(public_data_2),
  col = list(Cell_type = c("HSC"=my_palette1[2],"LMPP"=my_palette2[1],
                       "Cycling_LMPP"=my_palette2[2],"GMP1"=my_palette2[3],"GMP2"=my_palette2[4],
                       "Monocytes"=my_palette2[6],"DC"=my_palette2[7],
                       "CLP"=my_palette2[8], "PreB"=my_palette2[9],"MEP"=my_palette2[10],"Erythroid"=my_palette1[4],"Basophils"=my_palette2[5])),
  show_annotation_name = TRUE
  

)


pdf("Analysis/DEA/heatmap_public_data_sig_genes.pdf")
#cell_type > group > sample
Heatmap(public_data_2, name = "DEGs genes",
        col = structure(c(brewer.pal(n = 8, name = "YlGnBu")[c(7,5)], "gray", brewer.pal(n = 8, name = "YlOrRd")[c(4,7)]), names=c("-1","-0.5","0","0.5","1")),
        #col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 3), 
        cluster_rows = TRUE, cluster_columns = FALSE , use_raster=TRUE) 
dev.off()


## MERGE public data and our data in a heatmap

target_genes<-c("BRCA1", "BRCA2", "ATM", "ATR")



