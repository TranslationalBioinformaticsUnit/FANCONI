
########################################################

#### calculate binomial of each of the samples####

load("/datos_2/FANCONI/MERGE_4/heatmap_all.rda")

# select the data of each sample and split in two, one the fanconi sample and other the integrated data

fanconi <- new_dataset_2[,c(1:10)]
integrated <- new_dataset_2[,c(11:20)]

N <- nrow(fanconi_2004)
cells <- dim(fanconi_2004)[2]

a <- data.frame()
for (i in 1:N){
  for (j in 1:cells){
    a[i,j] <- ifelse( (fanconi[i,j]>0 & integrated[i,j]>0) | (fanconi[i,j]<0 & integrated[i,j]<0) | (fanconi[i,j]== integrated[i,j]),  1, 0 )
  }
}





Cell_type<-c("HSC","LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","MEP","Erythroid","Basophils")



pval_2004 <- c()
for (i in 1:ncol(a)){
     pval_2004[i] <- binom.test((table(a[,i])[2]), alternative="greater", N, p=1/2)$p.val
}

# adjust the p-value in each case

pval_adjusted_2004<-p.adjust(pval_2004, method = "bonferrani")

