### INTEGRATE ALL THE PUBLIC DATA (FA 1 AND 7 AND ALL THE HEALTHY ONES)
library(Seurat)
library(ggplot2)
library(ggpubr)

#>> Work directory and seed

setwd("/datos_2/FANCONI/Public_data")
set.seed(1234567)

FA_1 = Read10X(data.dir = "FA_1/outs/filtered_feature_bc_matrix/")
FA_1 <- CreateSeuratObject(counts = FA_1, project = "Fanconi", min.cells = 3, min.features = 400)
FA_7= Read10X(data.dir = "FA_7/FA_7_job/outs/filtered_feature_bc_matrix/")
FA_7 <- CreateSeuratObject(counts = FA_7, project = "Fanconi", min.cells = 3, min.features = 400)
Healthy_1=Read10X(data.dir = "Healthy_1/Healthy_1_job/outs/filtered_feature_bc_matrix/")
Healthy_1 <- CreateSeuratObject(counts = Healthy_1, project = "Fanconi", min.cells = 3, min.features = 400)
Healthy_2=Read10X(data.dir = "Healthy_2/Healthy_2_job/outs/filtered_feature_bc_matrix/")
Healthy_2 <- CreateSeuratObject(counts = Healthy_2, project = "Fanconi", min.cells = 3, min.features = 400)
Healthy_3=Read10X(data.dir = "Healthy_3/Healthy_3_job/outs/filtered_feature_bc_matrix/")
Healthy_3 <- CreateSeuratObject(counts = Healthy_3, project = "Fanconi", min.cells = 3, min.features = 400)
Healthy_4=Read10X(data.dir = "Healthy_4/Healthy_4_job/outs/filtered_feature_bc_matrix/")
Healthy_4 <- CreateSeuratObject(counts = Healthy_4, project = "Fanconi", min.cells = 3, min.features = 400)


#fanconi_corrected<-readRDS("/datos_2/FANCONI/DATA/Fanconi_CD34_MO.rds")
#fanconi_corrected$Patient<-"Fanconi_corrected"
#fanconi_corrected$Sample<-"Fanconi_corrected"
#fanconi_corrected$Project<-"Fanconi"
#
#
#healthy_our<-readRDS("/datos_2/FANCONI/DATA/seurat_young.rds")
#
#mo_193<-subset(healthy_our, subset=Patient=="mo193")
#mo_193$Sample<-"healthy_cima"
#mo_193$Project<-"Healthy"
#
#mo_194<-subset(healthy_our, subset=Patient=="mo194")
#mo_194$Sample<-"healthy_cima"
#mo_194$Project<-"Healthy"
#
#mo_196<-subset(healthy_our, subset=Patient=="mo196")
#mo_196$Sample<-"healthy_cima"
#mo_196$Project<-"Healthy"


##DEFINE PROJECT AND SAMPLE NAME TO THE ANALISYS

FA_1$Project<-"Fanconi"
FA_1$Sample<-"FA_1"
FA_1$Patient<-"FA_1"
FA_7$Project<-"Fanconi"
FA_7$Sample<-"FA_7"
FA_7$Patient<-"FA_7"
Healthy_1$Project<-"Healthy"
Healthy_1$Sample<-"Healthy_1"
Healthy_1$Patient<-"Healthy_1"
Healthy_2$Project<-"Healthy"
Healthy_2$Sample<-"Healthy_2"
Healthy_2$Patient<-"Healthy_2"
Healthy_3$Project<-"Healthy"
Healthy_3$Sample<-"Healthy_3"
Healthy_3$Patient<-"Healthy_3"
Healthy_4$Project<-"Healthy"
Healthy_4$Sample<-"Healthy_4"
Healthy_4$Patient<-"Healthy_4"


#>> INTEGRATION
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 10000 * 1024^2)

combined.list <- c(FA_1,FA_7,Healthy_1,Healthy_2,Healthy_3,Healthy_4)

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
                                     
saveRDS(fanconi.integrated, file="/datos_2/FANCONI/Public_data/integration_berria.RDS")

