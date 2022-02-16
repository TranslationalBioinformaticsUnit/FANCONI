## RUN AUCell for each fanconi patient

#----- LOAD libraries

library(AUCell)
library(GSEABase)
library(clusterProfiler)
library(ggplot2)


#------ Set work directory and create a folder for the analysis

setwd("/datos_2/FANCONI/AUCell_analysis/INTERSECTION_our_data")

#------ Load fanconi sample data and get expression matrix, normalized one


fanconi<-readRDS("/datos_2/FANCONI/DATA/fanconi_updated.RDS") # load the seurat object of our sample

exprMatrix<-as.matrix(fanconi@assays$RNA@data) # take the expression matrix, the normalized one

dim(exprMatrix)
#
fanconi@meta.data$CellType<- fanconi@meta.data$predicted.id


fanconi@meta.data$CellType<-factor(x=fanconi@meta.data$CellType, levels=c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils"))

fanconi<-SetIdent(fanconi, value="CellType")

levels(x = fanconi)<-c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils")
#
#
#
#
# CREATE CORRECTION GROUP
fanconi$FANCA<-fanconi@assays$RNA@data["FANCA",]
fanconi@meta.data$Correction <- as.character(fanconi@meta.data$FANCA > 0)


fanconi@meta.data$Correction[fanconi@meta.data$Correction==TRUE] <- "Corrected"
fanconi@meta.data$Correction[fanconi@meta.data$Correction==FALSE] <- "Uncorrected"


fanconi@meta.data$Correction<- as.factor(fanconi@meta.data$Correction)

##
##------ Load genesets using gmtfiles
#
#genes<-read.table("genes_pathway.txt")
##
#genes<-as.character(genes$x)
##
#gmtfile<-read.gmt("Cell_cycle.gmt.txt")
####
#genes<-gmtfile$gene
###
#geneSets<-GeneSet(genes, setName="Our_signature")

genes<-read.table("/datos_2/FANCONI/intersection_sig_genes_3_samples.txt", sep="\t") # load the list of genes to create a geneset
genes<-as.character(genes$x)
geneSets<-GeneSet(genes, setName="Intersection_our_4_samples") # create a geneset


## MYC pathway

#geneSets<-GeneSet(genes, setName="MYC_pathway")
#geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
#cbind(nGenes(geneSets))

#------ Calculat AUC score
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE) # generate a ranking, the genes are ordered alphabetically and the number is the position in the ranking

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings), aucMaxRank=2800) # calculate the AUC value, for default it only takes the top 5% genes, but we can modify this value with aucMaxRank
auc_per_cell <- getAUC(cells_AUC) # get the AUC value in a matrix for each one of the cells 
names(auc_per_cell_all) <- colnames(cells_AUC) # give the cell names
#------ Check the threshold and selected cells

rankings<-getRanking(cells_rankings) # get the matrix with all the rankings


#plot thresholds
set.seed(123)

cells_assignment <- AUCell_exploreThresholds(cells_AUC_all, plotHist=TRUE, assign=TRUE) # explore the automatic thresholds, not important in our case. We only check the scores in corrected and uncorrected cells

# check threshold and selected cells for fanconi pathway

cells_assignment$MYC_pathway$aucThr$thresholds # all the thresholds
#
cells_assignment$MYC_pathway$aucThr$selected # the selected one
#


# create a new variable with the score in the meta data

fanconi@meta.data$AUC_score<-rep(0, nrow(fanconi@meta.data))
fanconi@meta.data$AUC_score<-auc_per_cell_all[1,]



## boxplot per cell type

png("BM_healthy_score_correction.png")

ggplot(fanconi@meta.data, aes(x=CellType, y=AUC_score, fill=Correction)) + 
    geom_boxplot()+
    #scale_fill_manual(values=c("#D95F02"))+
    labs(title="Score of Cell cycle pathway cell type",x = "Cell_type", y = "score")
dev.off()
    
## boxplot per correction

png("BM_healthy_score_correction.png")

ggplot(prueba@meta.data, aes(x=CellType, y=Cluster1, fill=Correction)) + 
    geom_boxplot()+
    #scale_fill_manual(values=c("#D95F02"))+
    labs(title="Score of Cell cycle pathway by Correction level",x = "Cell_type", y = "score")
    
dev.off()
#
#metadata_BM_2002<-fanconi@meta.data
#metadata_BM_2002$Sample<-"BM_2002"
#
#save.image("BM_2002.RData")
#
##save.image("BM_healthy_mismatch_repair.RData")
###
#### MERGE all the metadata data frames and select the variable needed
##
#variables_to_keep<-c("CellType", "Correction", "AUC_score", "Sample")
###
#
##metadata_BM_2002<-metadata_BM_2002[,colnames(metadata_BM_2002)%in% variables_to_keep]
#metadata_BM_2004<-metadata_BM_2004[,colnames(metadata_BM_2004)%in% variables_to_keep]
#metadata_BM_2006<-metadata_BM_2006[,colnames(metadata_BM_2006)%in% variables_to_keep]
#metadata_BM_2008<-metadata_BM_2008[,colnames(metadata_BM_2008)%in% variables_to_keep]
#metadata_BM_healthy<-metadata_BM_healthy[,colnames(metadata_BM_healthy)%in% variables_to_keep]
##
#metadata_BM_all<-rbind(metadata_BM_2004, metadata_BM_2006, metadata_BM_2008, metadata_BM_healthy)
##metadata_BM_all$Tissue<-"BM"
###
###
###variables_pbmc<-c("Cell_type", "Correction", "AUC_score", "Sample")
###
###metadata_PBMC_2004<-metadata_PBMC_2004[,colnames(metadata_PBMC_2004)%in% variables_pbmc]
###metadata_PBMC_2006<-metadata_PBMC_2006[,colnames(metadata_PBMC_2006)%in% variables_pbmc]
###metadata_PBMC_2008<-metadata_PBMC_2008[,colnames(metadata_PBMC_2008)%in% variables_pbmc]
###metadata_PBMC_healthy<-metadata_PBMC_healthy[,colnames(metadata_PBMC_healthy)%in%  variables_pbmc]
####
###metadata_PBMC_all<-rbind(metadata_PBMC_2004, metadata_PBMC_2006, metadata_PBMC_2008, metadata_PBMC_healthy)
###metadata_PBMC_all$Tissue<-"PBMC"
####
###metadata_PBMC_all$CellType<-metadata_PBMC_all$Cell_type
####
###metadata_all$Sample = paste('BM', metadata_all$Sample, sep='_')
####colnames(metadata_all)<-c("CellType", "Correction", "AUC_score", "Sample", "Tissue", "Sample_2")
####
###metadata_all<-rbind(metadata_BM_all, metadata_PBMC_all)
####
###metadata_BM_all$Correction_tissue <- paste0(metadata_BM_all$Correction, "_", metadata_BM_all$Tissue)
###
###metadata_BM_all$Sample2<-metadata_BM_all$Sample
###metadata_BM_all$Sample2<-gsub("BM_", "", metadata_BM_all$Sample2)
###
###data_filtered$Sample<-gsub("PBMC_", "", data_filtered$Sample)
#####
####### plot the score by celltype
####metadata_all$Sample = paste('PBMC', metadata_all$Sample, sep="_")
#####
###    
#    p<-ggplot(data = metadata_BM_all, aes(x=Sample, y=AUC_score, fill=Correction)) + 
#         geom_boxplot() + facet_wrap(~CellType,ncol = 4)+
#    labs(title="Score of our intersection genes by cell type",x = "Cell_type", y = "score")
####
#pdf("BM_all_2002_score_sample.pdf", width=15, height=12)   
#   p + theme_bw()
#
#dev.off()
#####
####### plot the score by group by correction
####
####    
#   p<-ggplot(data = metadata_all, aes(x=Correction, y=AUC_score, fill=Correction)) + 
#        geom_boxplot() + facet_wrap(~Tissue,ncol = 4)+
#    labs(title="Score of our intersection genes by Correction level",x = "Cell_type", y = "score")
#  
#  pdf("ALL_TISSUES_TISSUE_boxplot_score_Correction.pdf", width=10)
####    p+theme_bw()
####dev.off()
####
####
###### plot density plot by sample
####
#p<-ggplot(prueba@meta.data, aes(x=Cluster1, fill=Correction)) +
# geom_density(alpha=0.4)+ facet_wrap(~Sample,ncol = 5)+
# labs(title="Score of our intersection genes")
# 
# p<-ggplot(metadata_BM_all, aes(x=AUC_score, fill="red")) +
# geom_density(alpha=0.4)+ facet_wrap(~Sample,ncol = 5)+
# labs(title="Score of our intersection genes")
# 
# 
####  
#  pdf("BM_ALL_2002_density_score.pdf", width=10)
#   p+theme_bw()
#dev.off()  
####  
#
#
p<-ggplot(fanconi@meta.data, aes(x=AUC_score, fill=Correction)) +
 geom_density(alpha=0.4)+ facet_wrap(~CellType,ncol = 4)+
 labs(title="Score of our intersection genes")
 
 
#   pdf("Density_celltype_healthy.pdf", width=10)
#   p+theme_bw()
#dev.off()  
#
#
#p1<- ggplot(metadata_BM_2006, aes(x=Sample, fill=Group)) + 
#  geom_bar(position="fill")

#-------- Check the effect of the number of genes selected to calculate the AUC, checking the number of genes that are also in our geneset 

summatory<-vector()
my_prop <- 0.2
for(i in 1:ncol(rankings)){
  prop <- nrow(rankings)*my_prop
  target_genes <- rownames(rankings)[rankings[,i]%in%c(1:prop)]
  a<-sum(target_genes %in% genes)
  summatory<-c(summatory,a)
}  
aa <- hist(summatory, breaks=50)

data_to_plot <- data.frame(mids=aa$mids,counts=aa$counts)
for(j in 1:(length(aa$breaks)-1)){  
  names <- colnames(rankings)[summatory>=aa$breaks[j] & summatory<=aa$breaks[j+1]]
  selected <- fanconi@meta.data[names,]
  data_to_plot$corrected[j] <- sum(selected$Correction=="Corrected")
  data_to_plot$uncorrected[j] <- sum(selected$Correction=="Uncorrected")
}


