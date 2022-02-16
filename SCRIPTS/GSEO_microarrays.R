######Compare microarray data with our data

library(ggplot2)
library(ComplexHeatmap)
library(xlsx)

##load our data information
setwd("/datos_2/FANCONI")

load("merge_samples_contrast_def.RDA")

## load microarray DEA data

setwd("/datos_2/FANCONI/GSE_microarray")

microarray<-read.xlsx("20210215_DEA_FANCA_vs_CONTROL_GSE16334.xlsx", sheetName="Hoja1")

sig_genes<-microarray[!microarray$sig_genes==0,]

sig_genes<-as.character(sig_genes$symbol)

write.table(sig_genes, "Mycroarray_dea_sig_genes.txt", sep="\t")

## create a data framee with our data and microarray data

fanconi_2006_sig<-fanconi_2006[sig_genes,]
uncorrected_2006_sig<-uncorrected_2006[sig_genes,]
corrected_2006_sig<-corrected_2006[sig_genes,]
data_2006<-cbind(fanconi_2006_sig, uncorrected_2006_sig, corrected_2006_sig)

fanconi_2004_sig<-fanconi_2004[sig_genes,]
uncorrected_2004_sig<-uncorrected_2004[sig_genes,]
corrected_2004_sig<-corrected_2004[sig_genes,]
data_2004<-cbind(fanconi_2004_sig, uncorrected_2004_sig, corrected_2004_sig)

fanconi_2008_sig<-fanconi_2008[sig_genes,]
uncorrected_2008_sig<-uncorrected_2008[sig_genes,]
corrected_2008_sig<-corrected_2008[sig_genes,]
data_2008<-cbind(fanconi_2008_sig, uncorrected_2008_sig, corrected_2008_sig)

data_fanconi<-cbind(data_2008,data_2006,data_2004)


colnames_to_grep<-c("genes_tendency", "symbol")

microarray_sig <- microarray[,colnames(microarray)%in%colnames_to_grep]

microarray_def<-microarray_sig[microarray$symbol%in%sig_genes,]

microarray_small<-microarray_def[microarray_def$genes_tendency==1|microarray_def$genes_tendency==-1,]


colnames(microarray_small)<-c("symbol", "microarray")
microarray_small$microarray<-microarray_small$microarray*(-1)

data_all<-merge(data_fanconi, microarray_small, by.x=0, by.y="symbol", all.y=TRUE, all.x=TRUE)



data_to_plot<-data_all[,2:98]

data_to_plot[is.na(data_to_plot)] <- 0

metadata_all <- data.frame(ID=colnames(data_to_plot),
                         cell_type=factor(c(rep(c("LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","MEP","Erythroid","Basophils"), times=3),rep(c("LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils"), times=3),rep(c("HSC","LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","MEP","Erythroid","Basophils"), times=3),"Microarray"),levels=c("HSC","LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils","Microarray")),
                          group=factor(c(rep(c("corrected_vs_uncorrected","healthy_vs_uncorrected","healthy_vs_corrected"), each=10),rep(c("corrected_vs_uncorrected","healthy_vs_uncorrected","healthy_vs_corrected"), each=11),rep(c("corrected_vs_uncorrected","healthy_vs_uncorrected","healthy_vs_corrected"), each=11),"microarray"), levels=c("corrected_vs_uncorrected","healthy_vs_uncorrected","healthy_vs_corrected","microarray")),
                          sample=factor(c(rep(c("sample_2008"),each=30),rep(c("sample_2006"),each=33),rep(c("sample_2004"),each=33), "microarray"),level=c("sample_2004","sample_2006","sample_2008","microarray")))


library(ComplexHeatmap)
library(RColorBrewer)



my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Cell_type = metadata_all$cell_type,
  Group=metadata_all$group,
  Sample=metadata_all$sample,
  col = list(Cell_type = c("HSC"=my_palette1[2],"LMPP"=my_palette2[1],
                       "Cycling_LMPP"=my_palette2[2],"GMP1"=my_palette2[3],"GMP2"=my_palette2[4],
                       "Monocytes"=my_palette2[6],"DC"=my_palette2[7],
                       "CLP"=my_palette2[8], "PreB"=my_palette2[9],"MEP"=my_palette2[10],"Erythroid"=my_palette1[4],"Basophils"=my_palette2[5], "Microarray"=my_palette1[6]),
              Group= c("corrected_vs_uncorrected"="#CD853F" ,"healthy_vs_uncorrected"= "#8B008B", "healthy_vs_corrected"="#6B8E23" , "microarray"=my_palette1[8]),
              Sample= c("sample_2008"="#FF6347", "sample_2006"="#6495ED", "sample_2004"="#008080", "microarray"=my_palette1[5])),
  show_annotation_name = TRUE
)
#plot
#jpeg("DEA/MM_C/heatmap_binary_sig_genes_MM_Control.jpg", width=5000, height=5000, res=600)

pdf("microarray_signature_our_data.pdf")
#cell_type > group > sample
Heatmap(data_to_plot, name = "DEGs genes",
        col = structure(c(brewer.pal(n = 8, name = "YlGnBu")[c(7,5)], "gray", brewer.pal(n = 8, name = "YlOrRd")[c(4,7)]), names=c("-1","-0.5","0","0.5","1")),
        #col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 3), 
        cluster_rows = TRUE, cluster_columns = FALSE, column_order=order(metadata_all$group,metadata_all$cell_type, metadata_all$sample) , use_raster=TRUE) 
dev.off()