##################################`
# 06/10/2020
#INTEGRATION HEALTHY AND FANCONI##
# MIREN LASAGA
##################################

library(Seurat)
library(ggplot2)
#library(ggpubr)

#>> Work directory and seed

setwd("/datos_2/FANCONI")
set.seed(1234567)

fanconi<-readRDS("/datos_2/FANCONI/DATA/Fanconi_2002/fanconi_names.RDS")

fanconi@meta.data$CellType<- fanconi@meta.data$predicted.id


fanconi@meta.data$CellType<-factor(x=fanconi@meta.data$CellType, levels=c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils"))

fanconi<-SetIdent(fanconi, value="CellType")

levels(x = fanconi)<-c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils")

fanconi$Sample<-"Fanconi_2002"


Healthy<-readRDS("DATA/mo268_names.rds")
Healthy$Sample<-"Healthy"




#>> INTEGRATION
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 10000 * 1024^2)

combined.list <- c(fanconi,Healthy)

for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], verbose = FALSE,vars.to.regress = c("percent.mt","nFeature_RNA"))
}

fanconi.features <- SelectIntegrationFeatures(object.list = combined.list, nfeatures = 3000)


combined.list <- PrepSCTIntegration(object.list = combined.list, anchor.features = fanconi.features, 
                                    verbose = FALSE)

combined.anchors <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT", 
                                           anchor.features = fanconi.features, verbose = FALSE)


fanconi.integrated <- IntegrateData(anchorset = combined.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)
                                     
saveRDS(fanconi.integrated, file="DATA/Fanconi_2002/integration.RDS")