##################################`
# 06/10/2020
#INTEGRATION HEALTHY AND FANCONI##
# MIREN LASAGA
##################################

library(Seurat)
library(ggplot2)
library(ggpubr)

#>> Work directory and seed

setwd("/datos_2/FANCONI")
set.seed(1234567)

#Fanconi<-readRDS("DATA/Fanconi_CD34_MO.rds")
#Fanconi@meta.data$CellType<- Fanconi@meta.data$predicted.id
#
#
#Fanconi@meta.data$CellType<-factor(x=Fanconi@meta.data$CellType, levels=c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils"))
#
#Fanconi<-SetIdent(Fanconi, value="CellType")
#
#levels(x = Fanconi)<-c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils")
#
#Fanconi@meta.data$Project<-"FANCONI"
#
#Fanconi@meta.data$Correction <- as.character(Fanconi@meta.data$FANCA > 0 )
#
#
#Fanconi@meta.data$Correction[Fanconi@meta.data$Correction==TRUE] <- "Corrected"
#Fanconi@meta.data$Correction[Fanconi@meta.data$Correction==FALSE] <- "Not_corrected"
#Fanconi@meta.data$Correction<- as.factor(Fanconi@meta.data$Correction)
#
#
#Healthy<-readRDS("DATA/seurat_young.rds")
#Healthy@meta.data$Project<-"HEALTHY"
#Healthy@meta.data$Correction<-"Healthy"

#>> SPLIT HEALTHY SAMPLE IN THE THREE INDIVIDUALS

#MO_193<-subset(Healthy, subset=Patient=="mo193")
#MO_194<-subset(Healthy, subset=Patient=="mo194")
#MO_196<-subset(Healthy, subset=Patient=="mo196")


sp_healthy<-readRDS("/datos_2/FANCONI/PUBLIC_SP/sp_public_updated.rds")
sp_healthy$Sample<-"Healthy"
sp_2006<-readRDS("/datos_2/FANCONI/Fanconi_SP/sp_2006.rds")
sp_2006$Sample<-"2006"

sp_2004<-readRDS("/datos_2/FANCONI/DATA/Fanconi_2004/SP/sp_2004_updated.rds")
sp_2004$Sample<-"2004"

sp_2008<-readRDS("/datos_2/FANCONI/DATA/Fanconi_2008/SP/sp_2008_updated.rds")
sp_2008$Sample<-"2008"


#>> INTEGRATION
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 10000 * 1024^2)

combined.list <- c(sp_healthy,sp_2006,sp_2004,sp_2008)

for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], verbose = FALSE)
}

fanconi.features <- SelectIntegrationFeatures(object.list = combined.list, nfeatures = 3000)


combined.list <- PrepSCTIntegration(object.list = combined.list, anchor.features = fanconi.features, 
                                    verbose = FALSE)

combined.anchors <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT", 
                                           anchor.features = fanconi.features, verbose = FALSE)


fanconi.integrated <- IntegrateData(anchorset = combined.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)
                                     
saveRDS(fanconi.integrated, file="DATA/SP_integrated.RDS")
