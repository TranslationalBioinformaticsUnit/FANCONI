#INPUT:
# List of target GO
# Output form GSEA --> GSEA set object

#-------------------------------
#gsea_object <- a
#path <- gsea_object@result$ID
path <- c("GO:0006298", "GO:0032201", "GO:0000723", "GO:0042769",  "GO:0032200", "GO:0032543", "GO:0000075", "GO:0006282", "GO:0009260","GO:0006302",
 "GO:0044839", "GO:0044786", "GO:0006275", "GO:0006301", "GO:0031570", "GO:0032781", "GO:0000724","GO:0006302", "GO:0006303", "GO:0007569","GO:0000086", "GO:0010971", "GO:0010972", "GO:0036297","GO:0052770")
 
# path <- c("GO:0006298", "GO:0032201", "GO:0000723", "GO:0042769",  "GO:0032200", "GO:0032543", "GO:0071166", "GO:0000075", "GO:0006282", "GO:0009260", "GO:0046034","GO:0006302", "GO:0033260",
# "GO:0044839", "GO:0044786", "GO:0042775", "GO:0050657", "GO:0000731",
#  "GO:0006275", "GO:0006301", "GO:0031570", "GO:0032205", "GO:0032781",
# "GO:0044819", "GO:0006977", "GO:0000724","GO:0000729","GO:0006302", "GO:0006303", "GO:0007569","GO:0000086", "GO:0010971", "GO:0010972", "GO:0036297","GO:0052770", "GO:0000302")
#---------------------------------



###### OUR_RIDGEPLOT FUNCTION ##############

our_ridgeplot = function(gsea_object, path, plot=TRUE){

if(!is.character(path)){
 stop("Error: vector with target pathways needed!")
}

if(!class(gsea_object)=="gseaResult"){
 stop("Error: gsea_object should be gseaResult class!")
}


require(ggridges)
library(BuenColors)

path_genes <- vector()
fc_path_genes <- vector()
path_id <- vector()
path_description <- vector()
pval <- vector() 

for(i in 1:length(path)){
    #select path position
    sel <- which(gsea_object@result$ID==path[i])
       
    #genes from pathway
    target_genes <- unlist(strsplit(gsea_object@result$core_enrichment[sel],"/"))
    path_genes <- c(path_genes,target_genes)
    
    #fch genes from pathway
    fc_target_genes <- gsea_object@geneList[target_genes]
    fc_path_genes <- c(fc_path_genes,fc_target_genes)
    
    #pathway ID
    path_id <- c(path_id,rep(path[i],length(target_genes)))
    
    #pathway description
    path_description <- c(path_description,rep(gsea_object@result$Description[sel],length(target_genes)))
    
    #pathway ajdusted p-valued
    pval <- c(pval,rep(gsea_object@result$p.adjust[sel],length(target_genes)))
    
    #control output
    cat(i,"-",path[i],": ",length(target_genes),"_",length(fc_target_genes),"\n")
}

#summary results table
data_for_ridgeplot <- data.frame(path_id,path_description,pval,path_genes,fc_path_genes)
data_for_ridgeplot$log10pval <- log10(data_for_ridgeplot$pval)*-1

#get_order
data_for_ridgeplot$path_description <- factor(data_for_ridgeplot$path_description, levels=gsea_object@result$Description[gsea_object@result$ID%in%path][order(gsea_object@result$p.adjust[gsea_object@result$ID%in%path], decreasing=TRUE)])

p<-"no_plot"
if(plot==TRUE){
p <- ggplot(data_for_ridgeplot, aes(x = fc_path_genes, y = path_description, fill=log10pval)) + 
      geom_density_ridges2() +
      scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits=c(0,3))
print(p)
}

#brewer_orange

return(list(results=data_for_ridgeplot, plot=p))
}

#----------------------------------------------------------------------------

#How to run:
#my_output <- our_ridgeplot(gsea_object, path, plot=TRUE)
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

  


###### OUR_RIDGE MULTIPLOT FUNCTION ##############
  
  
our_ridge_multiplot = function(list_gsea, path, gsea_names, order_by=NULL){

my_output_merged_list <- vector(mode="list", length=length(list_gsea))
my_output_merged_dataframe <- vector()

for(i in 1:length(list_gsea)){
  my_output <- our_ridgeplot(list_gsea[[i]], path, plot=FALSE)
  my_output$results$patient <- rep(gsea_names[i],nrow(my_output$results))
  my_output_merged_list[[i]] <- my_output$results
  my_output_merged_dataframe <- rbind(my_output_merged_dataframe, my_output$results)
}

if(!is.null(order_by)){
  sel <- which(gsea_names==order_by)
  my_output_merged_dataframe$path_description <- factor(my_output_merged_dataframe$path_description, levels=levels(my_output_merged_list[[sel]]$path_description))
}
  
p <- ggplot(my_output_merged_dataframe, aes(x = fc_path_genes, y = path_description, fill=log10pval)) + 
  geom_density_ridges(scale = 1) + facet_wrap(~factor(patient, levels=c("Corrected_uncorrected_2002", "Healthy_uncorrected_2002", "Corrected_uncorrected_2004","Healthy_uncorrected_2004","Corrected_uncorrected_2006","Healthy_uncorrected_2006","Corrected_uncorrected_2008", "Healthy_uncorrected_2008")))+
  scale_fill_gradientn(colors = jdb_palette("solar_rojos"), limits=c(0,3))+
  theme(text = element_text(size=15))

print(p)
return(list(my_output_merged_list,p))
}


#----------------------------------------------------------------------------

#How to run:
#my_output_2 <- our_ridge_multiplot(list_gsea, path, gsea_names, order_by=NULL)


  
  
  