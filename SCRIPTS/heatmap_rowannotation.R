library(ComplexHeatmap)
library(RColorBrewer)

setwd("C:/Users/transbio/Desktop/FANCONI/HEATMAPS")

load("heatmap_all_contrasts_def.RDA")

## get row names of target genes

target_TF <- c("DNMT1","MCM2","PCNA","MCM10","MCM7","UBE2T","DKC1","POLR3K","PARP1","FANCI","RAD51C","MSH6","RAD51AP1", "BRCA2", "CDK4", "MCM3", "MCM6",
               "RPA3", "HLA-A", "CDK1", "MCM5", "HLA-C", "HLA-E", "MCM4","UBE2C", "BRCA1")








target_TF_position <- vector()
for(i in 1:length(target_TF)){
  target_TF_position <- c(target_TF_position,grep(target_TF[i],rownames(new_dataset_2)))
}
rha = rowAnnotation(foo = anno_mark(at = target_TF_position,
                                    labels = rownames(new_dataset_2)[target_TF_position],
                                    labels_gp = gpar(fontsize = 6)))
## top annotation

metadata_all <- data.frame(ID=colnames(merge),
                           cell_type=factor(colnames(merge),levels=c("HSC","LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils")),
                           group=factor(c(rep(c("corrected_vs_uncorrected","healthy_vs_uncorrected","healthy_vs_corrected"), each=10),rep(c("corrected_vs_uncorrected","healthy_vs_uncorrected","healthy_vs_corrected"), each=11),rep(c("corrected_vs_uncorrected","healthy_vs_uncorrected","healthy_vs_corrected"), each=11)), levels=c("corrected_vs_uncorrected","healthy_vs_uncorrected","healthy_vs_corrected")),
                           sample=factor(c(rep(c("sample_2008"),each=30),rep(c("sample_2006"),each=33),rep(c("sample_2004"),each=33)),level=c("sample_2004","sample_2006","sample_2008")))


ha1 <- HeatmapAnnotation(
  
  Group=metadata_all$group,
  Cell_type = metadata_all$cell_type,
  Sample=metadata_all$sample,
  col = list(
             Group= c("corrected_vs_uncorrected"="#CD853F" ,"healthy_vs_uncorrected"= "#8B008B", "healthy_vs_corrected"="#6B8E23" ),
             Cell_type =c("HSC"=my_palette1[2],"LMPP"=my_palette2[1],
               "Cycling_LMPP"=my_palette2[2],"GMP1"=my_palette2[3],"GMP2"=my_palette2[4],
               "Monocytes"=my_palette2[6],"DC"=my_palette2[7],
               "CLP"=my_palette2[8], "PreB"=my_palette2[9],"MEP"=my_palette2[10],"Erythroid"=my_palette1[4],"Basophils"=my_palette2[5]),
             Sample= c("sample_2008"="#FF6347", "sample_2006"="#6495ED", "sample_2004"="#008080")),
  show_annotation_name = TRUE
)

merge_2<-merge[intersection,]


pdf("new_selection_heatmap.pdf")
#cell_type > group > sample
Heatmap(new_dataset_2, name = "DEGs genes",
        col = structure(c(brewer.pal(n = 8, name = "YlGnBu")[c(7,5)], "gray", brewer.pal(n = 8, name = "YlOrRd")[c(4,7)]), names=c("-1","-0.5","0","0.5","1")),
        #col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        right_annotation = rha,
        show_column_names = FALSE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 3), 
        cluster_rows = TRUE, cluster_columns = FALSE, column_order=order(metadata_all$group,metadata_all$cell_type, metadata_all$sample) , use_raster=TRUE) 
dev.off()


#cell_type > group > sample
Heatmap(merge_2, name = "DEGs genes",
        col = structure(c(brewer.pal(n = 8, name = "YlGnBu")[c(7,5)], "gray", brewer.pal(n = 8, name = "YlOrRd")[c(4,7)]), names=c("-1","-0.5","0","0.5","1")),
        #col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = TRUE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 3), 
        cluster_rows = TRUE, cluster_columns = FALSE, column_order=order(metadata_all$group, metadata_all$cell_type, metadata_all$sample) , use_raster=TRUE) 

