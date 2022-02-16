library(Seurat)
library(ggplot2)
#library(ggpubr)

#>> Work directory and seed

setwd("/datos_2/FANCONI")
set.seed(1234567)

#Fanconi_2008<-readRDS("DATA/Fanconi_2008/CD34/fanconi.RDS")
#
#Fanconi_2008$Sample<-"fanconi_2008"
#
#Fanconi_2004<-readRDS("DATA/fanconi_eltrombopag.RDS")
#Fanconi_2004$Sample<-"fanconi_2004"
#
#Fanconi_2006<-readRDS("DATA/fanconi_updated.RDS")
#
#Fanconi_2006$Sample<-"fanconi_2006"
#
#Fanconi_2002<-readRDS("DATA/fancon_2002.rds")
#
#Fanconi_2002$Sample<-"fanconi_2002"



Fanconi<-readRDS("DATA/4_integrated.RDS")
Fanconi$Project<-"Fanconi"

Healthy<-readRDS("DATA/mo268_names.rds")
Healthy$Project<-"Healthy"


#>> INTEGRATION
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 10000 * 1024^2)

combined.list <- c(Fanconi, Healthy)

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
                                     
saveRDS(fanconi.integrated, file="DATA/4_fanconi_healthy.RDS")