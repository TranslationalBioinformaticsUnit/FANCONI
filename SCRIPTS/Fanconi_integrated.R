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

## Number of corrected and uncorrected cells by cell type and sample

df<-read.csv("/datos_2/FANCONI/Fanconi_integrated/data_to_barplot.csv", sep=";")
df$Sample<-as.factor(df$Sample)

df$Cell_type<-factor(x=df$Cell_type, levels=c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils"))

q<-ggplot(df,aes(x=Sample,y=Count, fill=Correction))+
  geom_bar(stat = "identity",color="white")+
  facet_wrap(~Cell_type, nrow=1) +scale_y_break(c(18, 21)) 
  


scale_x_continuous(breaks = seq(-3, 4, 0.2)) 
pdf("Correction_by_celltype_def.pdf", width=13)
p1<-q + theme_bw()
p1+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) bp + coord_flip()

require(plotrix)
gap.barplot( as.matrix(df$Count), 
             beside = TRUE, col=as.factor(df$Correction),
             gap=c(1800,7000)), 
             ytics=c(0,3000,6000,9000,24000,25200,26400) )

dev.off()

## Compare DEGs of the 3 patients (use likert library)

load("merge_samples_contrast_def.RDA")

## in this RDA we have all the data for each one of the patients (fanconi_2004, fanconi_2006, fanconi_2008)


# create all the columns for each of the samples

fanconi_2006$HSC<-rep(0, nrow(fanconi_2006))
fanconi_2004$PreB<-rep(0, nrow(fanconi_2004))
fanconi_2008$HSC<-rep(0, nrow(fanconi_2008))
fanconi_2008$PreB<-rep(0, nrow(fanconi_2008))


## talke
sig_genes<-as.character(sig_genes$x)
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


## plot fanconi pathway for monocytes in each sample

setwd("Fanconi_integrated/KEEGplot/")

Sample<-c("2004","2006","2008") 



for (i in 1:length(Sample)){
  DE_genes<-read.table(paste0("Monocytes_",Sample[i],".txt"))
  png(paste0("Fanconi_integrated/Fanconi_KEGG_pathway_",Sample[i],".png"))
  pathview(gene.data = DE_genes, species = "hsa", pathway.id = "hsa03460", gene.idtype = "SYMBOL", out.suffix =Sample[i],kegg.native=FALSE)
 )
    dev.off()
  }



