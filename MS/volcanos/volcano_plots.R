library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)


KO1_LTED <- read.delim("KO-1G4_vs_LTED_significant_all.txt")
KO3_LTED <- read.delim("KO-3G4_vs_LTED_significant_all.txt")
WT_LTED <- read.delim("LTED_vs_WT_significant_all.txt")


volcano <- function(data){
  volcano_data <- data.frame(cbind(data[,23],data[,20],data[,18]))
  colnames(volcano_data) <- c("gene_name","log2fc","pvalue")
  volcano_data <- na.omit(volcano_data)
  volcano_data$log2fc <- as.numeric(volcano_data$log2fc)
  volcano_data$pvalue <- as.numeric(volcano_data$pvalue)
 
  volcano_data[1:length(volcano_data[,1]),4] <- "non-significant"
  indx_up <- which(((volcano_data$log2fc > 0.5) & (volcano_data$pvalue > 1.12)))
  volcano_data[indx_up,4] <- "up-regulated"
  
  indx_down <- which(((volcano_data$log2fc < -0.5) & (volcano_data$pvalue > 1.12)))
  volcano_data[indx_down,4] <- "down-regulated"
  colnames(volcano_data)[4] <- "regulation"
  return(volcano_data)
}


top_gene_names <- function(data,number_of_genes){
  data <- data[data$pvalue > 2 ,]
  sorted_data <- data %>% arrange(log2fc)
  head <- head(sorted_data,number_of_genes)
  sorted_data <- na.omit(sorted_data)
  tail <- tail(sorted_data,number_of_genes)
  top_genes <- rbind(head,tail) 
  return(top_genes)
} 

squash_axis <- function(from, to, factor) { 
  # A transformation function that squashes the range of [from, to] by factor on a given axis 
  
  trans <- function(x) { 
    
    # get indices for the relevant regions
    isq <- x > from & x < to & !is.na(x)
    ito <- x >= to & !is.na(x)
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor & !is.na(x)
    ito <- x >= from + (to - from)/factor & !is.na(x)
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(scales::trans_new("squash_axis", trans, inv, domain = c(from, to)))
}

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(family="Times New Roman"),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = "darkgray"),
            axis.title = element_text(face = "bold",size = 9.5),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.text.x = element_text(angle = 45, vjust = 0.5),
            axis.title.x = element_text(vjust = -0.2),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(),
            plot.margin=unit(c(3,3,3,3),"cm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

volcano_plot <- function(data_volcano,top_genes,title_n){
 
  p1 <- ggplot(data_volcano,aes(x=log2fc,y=pvalue)) + geom_point(aes(colour = regulation)) 
  p1 <- p1 + geom_text_repel(data = top_genes, aes(label = gene_name), 
                             alpha = 0.85,max.overlaps=12,min.segment.length = Inf,family="Times New Roman") 
  p1 <- p1 + ylim(0,8) + xlim(-5,5)
  p1 <- p1 + geom_vline(xintercept=c(-0.5, 0.5), col="darkgrey",linetype="dotted") 
  p1 <- p1 + geom_hline(yintercept=-2, col="darkgrey",linetype="dotted") 
  p1 <- p1 + scale_color_manual(name = "Regulation",values = c("down-regulated" = "steelblue", "up-regulated" = "tomato", "non-significant" = "grey")) 
  p1 <- p1 + labs(title = title_n,x = expression("t-test Difference"), y = expression(-log[10]~"(pvalue)"))
  p1 <- p1 + theme_Publication()
  return(p1)
}

KO1_LTED_volcano <- volcano(KO1_LTED)
top_genes_KO1_LTED <- top_gene_names(KO1_LTED_volcano,6)
fig1<- volcano_plot(KO1_LTED_volcano,top_genes_KO1_LTED,"LTED vs KO 1G4")

KO3_LTED_volcano <- volcano(KO3_LTED)
top_genes_KO3_LTED <- top_gene_names(KO3_LTED_volcano,6)
top_genes_KO3_LTED <- rbind(top_genes_KO3_LTED,KO3_LTED_volcano[KO3_LTED_volcano$gene_name =="GLYATL1",])
fig2 <- volcano_plot(KO3_LTED_volcano,top_genes_KO3_LTED,"LTED vs KO 3G4")

WT_LTED_volcano <- volcano(WT_LTED)
top_genes_WT_LTED <- top_gene_names(WT_LTED_volcano,6)
fig3 <- volcano_plot(WT_LTED_volcano,top_genes_WT_LTED,"LTED vs MCF7 WT")

library(ggpubr)
arrange <- ggarrange(fig1,fig2,fig3,ncol = 3, nrow = 1,common.legend=TRUE)
ggsave("arrangedplot_recent.png", arrange, width = 40, height = 9)
dev.off()


#MASS versus RPPA

rppa_proteins <- data.frame()
all_exp <- read.table("DIA_results.csv",sep=",",header=TRUE)
proteins <- c("SRC","MAPK1","NFKB1","ESR1","AKT1","RAF1","PIK3CA","STAT3")
for (protein in (1:length(proteins))){
  protein_exp <- all_exp[which(all_exp$Gene == proteins[protein]),c(4:18)]
  rppa_proteins <- rbind(rppa_proteins,cbind(proteins[protein],protein_exp))
}
rownames(rppa_proteins) <- rppa_proteins[,1]
rppa_proteins <- rppa_proteins[,-1]

sample_names <- data.frame(colnames(rppa_proteins))
colnames(sample_names) <- "sample"
sample_names <- sample_names %>% mutate(sample= gsub("_"," ",sample))
colnames(rppa_proteins) <- sample_names$sample

Heatmap(log(rppa_proteins), name = "log(MS Intensity)", rect_gp = gpar(col = "white", lwd = 1),
        column_title = "RPPA Expression of 11 Target Proteins in Knockdown Samples",
        cluster_columns =  F,
        column_names_rot = 45, 
        column_names_side = "top",
        column_dend_side = "bottom",
        col=brewer.pal(n = 6, name ="OrRd"),
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        heatmap_legend_param = list(direction = "horizontal"))


#### Enrichment #####3

KO1_LTED_volcano <- data.frame(volcano(KO1_LTED))
KO1_LTED_volcano <- na.omit(KO1_LTED_volcano)
KO1_LTED_volcano <- KO1_LTED_volcano[- grep(";", KO1_LTED_volcano$gene_name),]
write.table(data.frame(cbind(KO1_LTED_volcano$gene_name,KO1_LTED_volcano$log2fc)),"ms_enrich.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
KO1_LTED_volcano <- KO1_LTED_volcano %>% arrange(KO1_LTED_volcano,log2fc)

KO3_LTED_volcano <- volcano(KO3_LTED)
KO3_LTED_volcano <- na.omit(KO3_LTED_volcano)
KO3_LTED_volcano <- KO3_LTED_volcano[- grep(";", KO3_LTED_volcano$gene_name),]
KO3_LTED_volcano <- KO3_LTED_volcano %>% arrange(KO3_LTED_volcano,log2fc)
write.table(data.frame(cbind(KO3_LTED_volcano$gene_name,KO3_LTED_volcano$log2fc)),"ms_enrich_ko3_lted.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

WT_LTED_volcano <- volcano(WT_LTED)
WT_LTED_volcano <- na.omit(WT_LTED_volcano)
WT_LTED_volcano <- WT_LTED_volcano[- grep(";", WT_LTED_volcano$gene_name),]
WT_LTED_volcano <- WT_LTED_volcano %>% arrange(WT_LTED_volcano,log2fc)
write.table(data.frame(cbind(WT_LTED_volcano$gene_name,WT_LTED_volcano$log2fc)),"ms_enrich_wt_lted.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)


