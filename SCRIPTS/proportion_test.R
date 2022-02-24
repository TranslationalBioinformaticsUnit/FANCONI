## load fanconi data with the 4 samples integrated

fanconi<-readRDS("E:/FANCONI_analysis/DATA/integration_4_samples/4_integrated.RDS")


## create a data frame with total number of cells and corrected cells


fanconi$CellType<-factor(x=fanconi$CellType, levels=c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils"))
data<-data.frame(cell_type=rep(c(levels(as.factor(fanconi$CellType))),4))


data$Sample<-factor(c(rep(c("fanconi_2002"),each=12),rep(c("fanconi_2004"),each=12),rep(c("fanconi_2006"),each=12),rep(c("fanconi_2008"),each=12)),level=c("fanconi_2002","fanconi_2004","fanconi_2006","fanconi_2008"))
data$celltype_corrected<-rep(0, nrow(data))
data$total_corrected<-rep(0, nrow(data))
data$celltype_uncorrected<-rep(0, nrow(data))
data$total_uncorrected<-rep(0, nrow(data))
data$pvalue<-rep(0, nrow(data))


celltype<-as.character(data$cell_type)
sample<-as.character(data$Sample)

for(i in 1:length(celltype)){

data$celltype_corrected[i]<-sum(fanconi$CellType==celltype[i] & fanconi$Sample==sample[i] & fanconi$Correction=="Corrected")
data$total_corrected[i]<-sum(fanconi$Sample==sample[i] & fanconi$Correction=="Corrected")
data$celltype_uncorrected[i]<-sum(fanconi$CellType==celltype[i] & fanconi$Sample==sample[i] & fanconi$Correction=="Uncorrected")
data$total_uncorrected[i]<-sum(fanconi$Sample==sample[i] & fanconi$Correction=="Uncorrected")

res<-prop.test(x=c(data$celltype_uncorrected[i],data$celltype_corrected[i]),n=c(data$total_uncorrected[i],data$total_corrected[i]))
data$pvalue[i]<-res$p.value
}


data$p.adjust<-rep(0, nrow(data))

data$p.adjust<-p.adjust(data$pvalue, method="bonferroni")



data$corrected[i]<-sum(fanconi$CellType==celltype[i] & fanconi$Correction=="Corrected")