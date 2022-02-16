#################### CLUSTERS OF PATHWAYS OBTAINED FROM ENRICHMENT ANALYSIS ####################

library(GeneSetCluster)   # https://github.com/TranslationalBioinformaticsUnit/GeneSetCluster/wiki




# ---------------------

# STEP 1.1: LOAD THE DATA

setwd("/Users/aurdangarin/Documents/STATegRa GSEA/Results/Results_NPC_all_GBM")

pathways_all <- read.table("Pathways_rna_newrna.csv", sep=";", header=TRUE)
pathways_rna <- read.table("Pathways_rna.csv", sep=";", header=TRUE)
pathways_new.rna <- read.table("Pathways_new_rna.csv", sep=";", header=TRUE)




# STEP 1.2: MERGE THE DATA

pathways <- pathways_all

pathways$rna <- ifelse(pathways$ID %in% pathways_rna$ID, "yes", "no")
pathways$new.rna <- ifelse(pathways$ID %in% pathways_new.rna$ID, "yes", "no")
pathways$all <- ifelse(pathways$ID %in% pathways_all$ID, "yes", "no")

head(pathways)





# ---------------------

# STEP 2: CREATE PATHWAY OBJECT

object <- ObjectCreator(Pathways = pathways$ID, Molecules = pathways$geneID, Groups=rep("union", nrow(pathways)),
                        structure = "Entrez", Type="", sep="/", Source="enricher", organism="org.Hs.eg.db")

# Groups, structure, Type, Source, organism -> strings with useful information for us


ShowExperimentdata(Object = object)
ShowMeta(Object = object)





# ---------------------

# STEP 3: MAKE CLUSTERS (compute RR)

combine <- CombineGeneSets(Object = object, display="Condensed")
combine@Data.RR
combine@metadata
combine@PData
combine@plot




# STEP 3.1: SELECT THE OPTIMAL NUMBER OF CLUSTERS

library(ggplot2)
OptimalGeneSets(object = combine, method = "silhouette", max_cluster= 20, cluster_method = "kmeans", main= "Kmeans for 20 clusters")




# STEP 3.2: CLUSTER GENE SETS

clusters <- ClusterGeneSets(Object = combine, 
                            clusters = 5, 
                            method = "kmeans", 
                            order = "Sum.RR", 
                            molecular.signature = "All")

clusters@Data.RR
clusters@Data
clusters@metadata
clusters@plot


# STEP 3.3: PLOT THE RR 

# Prepare the annotations data

annot <- pathways[, colnames(pathways) %in% c("ID", "rna", "new.rna", "all") ]

annot$ID <- as.character(annot$ID)
annot$ID <- paste0("union_", annot$ID)

annot <- annot[match(rownames(clusters@Data.RR), annot$ID),]

annot$cluster <- as.factor(clusters@Data[[1]]$cluster)
rownames(annot) <- annot$ID
annot$ID <- NULL




# Plot RR

library(pheatmap)
library(RColorBrewer)

my_palette1 <- c(brewer.pal(11, "Set3")[c(1:12)])
my_palette2 <- c(brewer.pal(11, "PRGn")[c(1:11)])

ann_colors = list(
  all = c(yes="#80B1D3", no="#D9D9D9"),
  rna= c(yes="#80B1D3", no="#D9D9D9"),
  new.rna=c(yes="#80B1D3", no="#D9D9D9"), 
  cluster=c('1'="#8DD3C7", '2'="#FFFFB3", '3'="#BEBADA", '4'="#FB8072", '5'="#1B7837", '6'="#FDB462", '7'="#B3DE69", '8'="#FCCDE5", '9'="#A6DBA0", 
            '10'="#5AAE61", '11'="#CCEBC5", '12'="#A6CEE3"))



PlotGeneSets(Object = clusters, fontsize = 6,
             legend = T,
             annotation.mol=F, RR.max = 60,
             main="", annot.col.rows = annot, colors=ann_colors)




# ---------------------


# SAVE CLUSTERS AND PATHWAYS IN A .csv FILE

cl <- as.data.frame(clusters@Data)
head(cl)

cl <- cl[, c("Pathways", "Molecules", "cluster")]
colnames(cl) <- c("Pathways", "Genes", "Cluster")


colnames(pathways)



pathways <- pathways[, c("ID", "Description", "rna", "new.rna", "all")]
head(pathways)

colnames(pathways) <- c("Pathways", "Description", "rna", "new.rna", "rna.and.newrna")




pathways_cluster <- merge(pathways, cl, by="Pathways")
head(pathways_cluster)

pathways_cluster2 <- pathways_cluster[order(pathways_cluster$Cluster, decreasing=FALSE), ]
head(pathways_cluster2)



write.table(pathways_cluster2, file="/Users/aurdangarin/Documents/STATegRa GSEA/Results/Results_NPC_all_GBM/Clusters_pathways_all.csv", 
            sep=";", row.names = FALSE)
