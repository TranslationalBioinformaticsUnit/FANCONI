
########################################################

#### 
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


a <- data.frame()
for (i in 1:N){
  for (j in 1:cells){
    a[i,j] <- ifelse( fanconi_2004[i,j]>0 ,  1, 0)
  }
}
b<-vector()
for (i in 1:11){
    c<-
    b[i] <- cor.test(as.vector(a[,i]),as.vector(a[,i+11]), method="spearman", exact=FALSE)$p.val
  }


Cell_type<-c("HSC","LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","MEP","Erythroid","Basophils")



pval_2004 <- c()
for (i in 1:ncol(a)){
     pval_2004[i] <- binom.test((table(a[,i])[2]), alternative="greater", N, p=1/2)$p.val
}


