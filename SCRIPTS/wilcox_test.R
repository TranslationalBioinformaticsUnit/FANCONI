




## load fanconi data of the 4 samples integrated

fanconi<-readRDS("E:/FANCONI_analysis/DATA/integration_4_samples/4_fanconi_healthy.RDS")

cell_type<-c("HSC", "LMPP", "Cycling_LMPP", "GMP1", "GMP2", "Monocytes", "DC", "CLP", "PreB", "MEP", "Erythroid", "Basophils")
pvalue<-vector()

for(i in 1:length(cell_type)){

data<-subset(fanconi, subset=CellType==cell_type[i])

res<-wilcox.test(data$FANCA ~ data$Project, data=data@meta.data, alternative="two.sided")
#res<-aov (data$FANCA ~ data$Sample, data=data@meta.data)

pvalue<-c(pvalue, res$p.value)}

