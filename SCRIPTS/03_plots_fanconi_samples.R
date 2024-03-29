# load required libraries
library(Seurat)
library(ggplot2)
library(ggpubr)
library(likert)
library(grid)
library(RColorBrewer)

set.seed(1234567)
##set work directory
setwd("/datos_2/FANCONI")
#####
fanconi<-readRDS("DATA/fanconis_integrated.RDS")

fanconi@meta.data$CellType<-factor(x=fanconi@meta.data$CellType, levels=c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils"))

##FANCA EXPRESSION


corrected<-subset(fanconi, subset=FANCA>0)

p<-ggplot(corrected@meta.data, aes(x=CellType, y=FANCA, fill=Sample)) + 
    geom_boxplot()+
    scale_fill_manual(values=c("lightpink", "#008080","#6495ED", "#FF6347"))+
    labs(title="Expression of FANCA by cell type",x = "Cell_type", y = "FANCA")

pdf("FANCA_expression_fanconis.pdf", width=10)
p +theme_bw()
dev.off()

#read data
data<-read.csv("data_to_barplot.csv", sep=";")

#reorder the data

data$Cell_type<-factor(data$Cell_type, levels<-c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils")
)

#raw plot
p0 <- ggplot(data,aes(x=as.factor(Sample),y=Count, fill=Correction))+
  geom_bar(stat = "identity",color="white")+
  facet_wrap(~Cell_type, nrow=1) +
  ylim(0,9000) +
  theme_bw() 

p0 + scale_y_break(c(520,700), scales=0.10, space=0.01) + scale_y_break(c(1800,7400), scales=0.05, space=0.01)

pdf("barplot_number_cells_breaks.pdf", width=15)
p0 + scale_y_break(c(520,700), scales=0.10, space=0.01) + scale_y_break(c(1800,7400), scales=0.05, space=0.01)
dev.off()

## Compare DEGs of the 4 patients (use likert library)

load("/datos_2/FANCONI/MERGE_4/heatmap_all.rda")

## in this RDA we have all the data for each one of the patients (fanconi_2002, fanconi_2004, fanconi_2006, fanconi_2008)


# create all the columns for each of the samples

fanconi_2006$HSC<-rep(0, nrow(fanconi_2006))
fanconi_2004$PreB<-rep(0, nrow(fanconi_2004))
fanconi_2008$HSC<-rep(0, nrow(fanconi_2008))
fanconi_2008$PreB<-rep(0, nrow(fanconi_2008))


## talke
sig_genes<-rownames(union) # take the genes that are significant in at least one sample and one celltype
fanconi_2002_sig<-data_fanconi[sig_genes,]
fanconi_2004_sig<-fanconi_2004[sig_genes,]
fanconi_2006_sig<-fanconi_2006[sig_genes,]
fanconi_2008_sig<-fanconi_2008[sig_genes,]

general<-rbind(fanconi_2002_sig,fanconi_2004_sig, fanconi_2006_sig, fanconi_2008_sig)


group<-c(rep("d_fanconi_2002",nrow(fanconi_2002_sig)),rep("c_fanconi_2004",nrow(fanconi_2004_sig)),rep("b_fanconi_2006",nrow(fanconi_2006_sig)), rep("a_fanconi_2008",nrow(fanconi_2008_sig)))#,levels=c("fanconi_2002", "fanconi_2004", "fanconi_2006", "fanconi_2008")

levels(group)=c("fanconi_2008", "fanconi_2006", "fanconi_2004", "fanconi_2002")

for(i in 1:ncol(general)){
 general[,i] <- factor(general[,i],levels=c(-1,-0.5,0.5,1))
}

colnames(general)<-c("a_HSC", "b_LMPP","c_Cycling_LMPP","d_GMP1","e_GMP2","f_Monocytes","g_DC","h_CLP","i_PreB","j_MEP","k_Erythroid","l_Basophils")
celltype<-factor(x=colnames(general), levels=c("a_HSC", "b_LMPP","c_Cycling_LMPP","d_GMP1","e_GMP2","f_Monocytes","g_DC","h_CLP","i_PreB","j_MEP","k_Erythroid","l_Basophils"))
p <- our.likert(general, grouping=group)


source("/datos_2/FANCONI/SCRIPTS/our.likert.R")
pdf("likert_plot_UNION_genes.pdf", height=10, width=10)
our.likert.bar.plot(p,include.histogram = TRUE, color=c(brewer.pal(n = 8, name = "YlGnBu")[c(7,5)],brewer.pal(n = 8, name = "YlOrRd")[c(4,7)]))
dev.off()




