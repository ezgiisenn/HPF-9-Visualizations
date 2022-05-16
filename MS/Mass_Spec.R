library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)

#Mass Spec
dia_results <- read.table("DIA_results.csv",sep=",",header=TRUE,stringsAsFactors = F)

#quality control

#1. intensity distribution of the data
pdf(file="/Users/ezgisen/Documents/HPF9/Mass_Spec/plots_Ms_script/plots.pdf",paper="USr")

boxplot(dia_results[,c(4:19)],las=2,main="Boxplot of the MS DATA")
plot(density(dia_results$WT_r3,na.rm=T),xlim=c(0,5e6))

#log2 transform the data

dia_results[,4:19] <- log2(dia_results[,4:19])
samples <- data.frame("name"=colnames(dia_results[,4:19]))

samples <- samples  %>% mutate(name= gsub("_"," ",name))
colnames(dia_results)[c(4:19)] <- c(samples$name)
dia <- dia_results[,c(4:19)]

par(mar=(c(10,5,5,5)))
boxplot(dia,col = rep(c("lightblue","black","seagreen","tan1"),each=4),las=2,main="Boxplot of the MS DATA")

#protein IDs

par(mar=(c(10,5,5,5)))
protein_IDs <- apply(dia,2, function(x) {sum(!is.na(x))} ) 
barplot(protein_IDs,las=2,main="Number of Protein Groups Identified",
        ylab="number of protein",
        col= rep(c("lightblue","black","seagreen","tan1"),each=4))

#exclude WT_4
dia_results <- dia_results[,-19]

#how many values are there in each row for  each peptide group
dia_results$no <- apply(dia_results[,4:18],1, function(x) {sum(!is.na(x))} ) 
table(dia_results$no)/sum(table(dia_results$no))

#Scree Plot PCA 
par(mar=(c(10,5,5,5)))
dia_pca <- prcomp(t(dia_results[dia_results$no == 15,c(4:18)]),scale=TRUE)
var_exp <- dia_pca$sdev
var_expl <-((var_exp)*2/sum(var_exp)*2) * 100
barplot(var_expl,beside=TRUE,names.arg = colnames(dia_pca$rotation),
        main="Scree Plot",col="steelblue",ylab="(%)")

#PCA
par(mar=(c(10,10,5,5)))
plot(dia_pca$x[,1], dia_pca$x[,2],pc=16,
     col=rep(c("lightblue","black","seagreen","tan1"),each=4),
     ylab="PC2 (58.33 %)",xlab="PC1 (82.79 %)",main="PCA")

legend("topright",col= c("lightblue","black","seagreen","tan1"),pch=16,
     legend=c("KO 1G4","KO 3G4","LTED","WT"),bty="n")

#heatmap
par(mar=(c(10,5,5,5)))
sample_col <- data.frame("sample" = c(rep(c("KO_1G4","KO_3G4","LTED"),each =4),rep("WT",3)))
rownames(sample_col) = colnames(dia_results)[c(4:18)]
my_colour = list(
  sample = c(KO_1G4 = "lightblue",KO_3G4 = "black",LTED = "seagreen",WT = "tan1"))
pheatmap(dia_results[dia_results$no == 15,4:18],scale="row",show_rownames=FALSE,
         annotation = sample_col,
         annotation_colors = my_colour,
         cutree_cols= 4)
dev.off()

write.table(dia_results,"dia_results_filtered.csv", row.names=FALSE, sep=",")















