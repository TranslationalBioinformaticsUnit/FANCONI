#Barplot fanconi

library(ggplot2)
library(ggbreak)

#read data
data<-read.csv("data_to_barplot.csv", sep=";")

#reorder the data

data$Cell_type<-factor(data$Cell_type, levels<-c("HSC", "LMPP","Cycling_LMPP","GMP1","GMP2","Monocytes","DC","CLP","PreB","MEP","Erythroid","Basophils")
)

#raw plot
p0 <- ggplot(data,aes(x=as.factor(Sample),y=Count, fill=Correction))+
  geom_bar(stat = "identity",color="white")+
  facet_wrap(~Cell_type, nrow=1) +
  ylim(0,9000) +
  theme_bw() 


data_corr <- data[data$Correction=="Corrected",]


ggplot(data,aes(x=as.factor(Sample),y=rep(1,nrow(data)), size=Count, color=Count))+
  geom_point(alpha=0.5)+
  facet_wrap(~Cell_type, nrow=1) +
  theme_bw() 

require("BuenColors")

pdf("buble_number_cells.pdf", width=15)
ggplot(data,aes(x=as.factor(Sample),y=rep(1,nrow(data)), fill=Count))+
  geom_point(alpha=1, shape=22, color="white", size=8)+
  scale_fill_gradientn(colors = jdb_palette("brewer_spectra"))+
  facet_wrap(~Cell_type, nrow=1) +
  theme_bw() 
dev.off()
p0 + scale_y_break(c(520,700), scales=0.10, space=0.01) + scale_y_break(c(1800,7400), scales=0.05, space=0.01)
  
pdf("barplot_number_cells_breaks.pdf", width=15)
p0 + scale_y_break(c(520,700), scales=0.10, space=0.01) + scale_y_break(c(1800,7400), scales=0.05, space=0.01)
dev.off()

# plot percentages
p0 <- ggplot(data,aes(x=as.factor(Sample), y=Count, fill=Correction))+
  geom_bar(position="fill", stat="identity")+
  facet_wrap(~Cell_type, nrow=1) +
  theme_bw() 

## plot corrected and uncorrected separately and log transformed
p1 <- ggplot(data,aes(x=as.factor(Sample),y=Count, fill=Correction))+
  geom_bar(stat = "identity",color="white", position="dodge")+
  facet_wrap(~Cell_type, nrow=1) +
  ylim(0,9000) +
  theme_bw() + scale_y_continuous(trans="log10")

pdf("barplot_log10_dodge.pdf", width = 15)
p1
dev.off()

#plot with log transformation
#prepare data to plot
data_to_plot <- data.frame(cell_type=rep(unique(data$Cell_type),4), 
                            Sample=rep(unique(data$Sample), each=length(unique(data$Cell_type))))

for(i in 1:nrow(data_to_plot)){
  data_to_plot$count_corrected[i] <- data$Count[data$Cell_type==data_to_plot$cell_type[i] & data$Sample==data_to_plot$Sample[i] & data$Correction=="Corrected"]
  data_to_plot$count_uncorrected[i] <- data$Count[data$Cell_type==data_to_plot$cell_type[i] & data$Sample==data_to_plot$Sample[i] & data$Correction=="Uncorrected"]
}


data_to_plot$total <- data_to_plot$count_corrected + data_to_plot$count_uncorrected
data_to_plot$log10_total <- log10(data_to_plot$total)
data_to_plot$log10_count_corrected <- log10(data_to_plot$count_corrected)
data_to_plot$log10_count_uncorrected <- log10(data_to_plot$count_uncorrected)
data_to_plot$diff <- data_to_plot$log10_total-data_to_plot$log10_count_uncorrected


final_data_to_plot <- data.frame(cell_type=c(data_to_plot$cell_type,data_to_plot$cell_type),
                                 sample=c(data_to_plot$Sample,data_to_plot$Sample),
                                 correction=c(rep("uncorrected",nrow(data_to_plot)),rep("corrected",nrow(data_to_plot))),
                                 log10_value=c(data_to_plot$log10_count_uncorrected,data_to_plot$diff)
                                 )

final_data_to_plot$cell_type <- factor(final_data_to_plot$cell_type, levels = unique(final_data_to_plot$cell_type))

#Plot
p <- ggplot(final_data_to_plot,aes(x=as.factor(sample),y=log10_value, fill=correction))+
  geom_bar(stat = "identity",color="white")+
  facet_wrap(~cell_type, nrow=1)

p + scale_y_continuous(labels = 10^c(0:4))  
