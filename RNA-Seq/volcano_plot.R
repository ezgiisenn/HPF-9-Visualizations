library(readxl)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(scales)

data <- read_excel("20220511_wt_vs_LTED_DESeq2.xlsx")

for (gene_anno in (1:length(data$hgnc_symbol))){
  if (is.na(data$hgnc_symbol[gene_anno] == TRUE)){
    data$hgnc_symbol[gene_anno] <- data$ensembl_gene_id[gene_anno]
  }
}

volcano <- data.frame(cbind(data$hgnc_symbol,data$log2FoldChange,data$padj))


#make sure the values read from excel are in numeric or double
colnames(volcano) <- c("gene_name","log2fc","pvalue")
volcano[,2] <- as.numeric(volcano[,2])
volcano[,3] <- as.numeric(volcano[,3])

#remove na values to not to get a warning later
volcano <- volcano[!is.na(volcano$pvalue),]
volcano <- volcano[!is.na(volcano$log2fc),]

#annotate each gene as non-significant, up & down for further coloring
volcano[,"regulation"] <- "non-significant"
indx_up <- which(((volcano$log2fc > 0.5) & (volcano$pvalue < 0.01)))
volcano[indx_up,"regulation"] <- "up"

indx_down <- which(((volcano$log2fc < -0.5) & (volcano$pvalue < 0.01)))
volcano[indx_down,"regulation"] <- "down"
volcano <- volcano %>% arrange((log2fc), pvalue)

#see top upregulated genes
head(volcano) 

#see top downregulated genes
tail(volcano)

head <- volcano %>% filter(row_number() < 6)
tail <- volcano %>% filter(row_number() > length(volcano[,1])-5)
top_genes <- rbind(head,tail)


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
  return(trans_new("squash_axis", trans, inv))
}


pdf("volcano_plot_squashed.pdf",paper="USr")
p <- ggplot(volcano,aes(x=log2fc,y=-log10(pvalue))) + geom_point(aes(colour = regulation)) 
p <- p + geom_text_repel(data = top_genes, aes(label = gene_name), 
                         alpha = 0.85,max.overlaps=12,min.segment.length = Inf)
p <- p + coord_trans(y = squash_axis(30, 125, 8))
p <- p + theme_bw() 
p <- p + geom_vline(xintercept=c(-0.5, 0.5), col="darkgrey",linetype="dotted") + geom_hline(yintercept=-log10(0.01), col="darkgrey",linetype="dotted")
p <- p + scale_color_manual(name = "Regulation",values = c("down" = "royalblue1", "up" = "salmon", "non-significant" = "grey"))
p <- p + xlab(expression(log[2]~"Fold Change")) + ylab(expression(-log[10]~"(pvalue)"))
p <- p + ggtitle("Volcano Plot of DE Genes LTED vs. MCF7")
p <- p + theme(legend.position="bottom", legend.box = "horizontal", legend.title = element_blank(),
               plot.title = element_text(hjust=0.5,face = "bold"))
p
dev.off()

#heatmap
data_transformed <- data
data_transformed[,4:9] <- log2(data[,c(4,5,6,7,8,9)]+1)
data_transformed <- data_transformed %>% arrange((log2FoldChange), padj)
colnames(data_transformed)[4:9] <- c("MCF7 WT 2","MCF7 WT 3","MCF7 WT 4","MCF7 LTED 2",
                                     "MCF7 LTED 3","MCF7 LTED 4")

data_transformed_up <- data_transformed %>% filter(row_number() <= 20)
data_transformed_down <- data_transformed %>% filter(row_number() >= 24755)

data_heatmap_end <- rbind(data_transformed_up,data_transformed_down)
annotation <- data_heatmap_end$hgnc_symbol

data_end <- as.data.frame(data_heatmap_end[,4:9])
rownames(data_end) <- annotation
col_info <- data.frame(sample = rep(c("WT", "LTED"), c(3,3)))
rownames(col_info) <- colnames(data_end)

color <- list(sample = c(WT = "tan1", LTED="seagreen"))
pdf("heatmap1.pdf",paper="USr")
pheatmap(data_end, annotation_col = col_info,cluster_rows=FALSE, 
         cutree_rows = 2, 
         cutree_cols = 2,
         annotation_colors = color[1],
         gaps_row=20,
         main = "Top 40 Regulated Genes MCF7 WT vs MCF7 LTED")
dev.off()
annotation_colors = annotation_col=c('salmon','royalblue')
