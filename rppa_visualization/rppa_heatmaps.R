library(readxl)
library(pheatmap)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(grDevices)


data <- read_excel("z-scaled.xlsx")
data <- data.frame(data[,-1])
data <- data %>% mutate(sample_name= gsub("_"," ",sample_name))
rownames(data) <- data[,1]
data <- data[,-1]

suppressPackageStartupMessages(library(ComplexHeatmap))
ht_opt(heatmap_column_names_gp = gpar(fontfamily = "Times"), 
       heatmap_row_names_gp = gpar(fontfamily = "Times"), 
       heatmap_column_title_gp = gpar(fontfamily = "Times"),
       heatmap_row_title_gp = gpar(fontfamily = "Times"),
       legend_title_gp = gpar(fontfamily = "Times"),
       legend_labels_gp = gpar(fontfamily = "Times"))

Heatmap(t(data), name = "z-scaled log \nExpression", rect_gp = gpar(col = "white", lwd = 1),
        #column_title = "RPPA Expression of All Target Proteins in RPPA across Knockdown Samples",
        cluster_columns = FALSE, cluster_rows=FALSE,
        column_names_rot = 45, column_names_side = "top",
        #column_dend_side = "bottom",
        col=brewer.pal(n = 6, name ="OrRd"),
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        heatmap_legend_param = list(direction = "horizontal"))

draw(ht_list, heatmap_legend_side = "bottom")

data <- data.frame(data[-c(1,4,7,10),])
data <- data.frame(data[,-c(11,12)])

#ES medium state based Clustering
data <- data[c(1,2,5,6,3,4,7,8),]

#KD state based Clustering

pdf("rppa_heatmap_without_sample1",paper="USr")
par(mar = c(10, 10, 10, 10))
pheatmap(t(data))
dev.off()
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(family="Times New Roman"),
            panel.background = element_rect(colour = NA),
            plot.background = element_blank(),
            panel.border = element_blank(),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.text.x = element_text(angle = 45,hjust=1,size=16),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.line = element_line(colour="darkgray"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.position="none",
            plot.margin=unit(c(1,1,1,1),"cm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  

}

squash_axis <- function(from, to, factor) { 
  # A transformation function that squashes the range of [from, to] by factor on a given axis 
  
  # Args:
  #   from: left end of the axis
  #   to: right end of the axis
  #   factor: the compression factor of the range [from, to]
  #
  # Returns:
  #   A transformation called "squash_axis", which is capsulated by trans_new() function
  
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

##ESR Boxplot##
data <- read_excel("z-scaled.xlsx")
ESR <- data.frame(data$ESR)
rownames(ESR) <- rownames(data)
ESR <- cbind(ESR,data.frame(c(rep(c("KD ctrl ES-","KD ctrl ES +","KD ES-","KD ES+"),each=3))))
ESR <- cbind(ESR,data.frame(c("ES-","ES-","ES-","ES+","ES+","ES+",
                              "ES-","ES-","ES-","ES+","ES+","ES+")))
colnames(ESR)[c(1,2,3)] <- c("expression","sample","medium")

ggplot(ESR, aes(x=sample, y=expression,color=medium)) + 
  geom_boxplot() +
  stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
  theme_Publication() + 
  scale_color_manual(values=c("ES-" ="#F7C1BB", "ES+"="#F15152"))+
  labs(x="",y="z-scaled expression",title="ESR1")

#pAkt1 Boxplot####
pAkt1 <- data.frame(data$pAkt1)

rownames(pAkt1) <- rownames(pAkt1)
pAkt1 <- cbind(pAkt1,data.frame(c(rep(c("KD ctrl ES-","KD ctrl ES +","KD ES-","KD ES+"),each=3))))
pAkt1 <- cbind(pAkt1,data.frame(c("ES-","ES-","ES-","ES+","ES+","ES+",
                              "ES-","ES-","ES-","ES+","ES+","ES+")))
colnames(pAkt1)[c(1,2,3)] <- c("expression","sample","medium")

ggplot(pAkt1, aes(x=sample, y=expression,color=medium)) + 
  geom_boxplot() +
  stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
  theme_Publication() + 
  scale_color_manual(values=c("ES-" ="#F7C1BB", "ES+"="#F15152"))+
  labs(x="",y="z-scaled expression",title="pAkt1")

#EFGR Boxplot###

pAkt1 <- data.frame(data$STAT3)

rownames(pAkt1) <- rownames(pAkt1)
pAkt1 <- cbind(pAkt1,data.frame(c(rep(c("KD ctrl ES-","KD ctrl ES +","KD ES-","KD ES+"),each=3))))
pAkt1 <- cbind(pAkt1,data.frame(c("ES-","ES-","ES-","ES+","ES+","ES+",
                                  "ES-","ES-","ES-","ES+","ES+","ES+")))
colnames(pAkt1)[c(1,2,3)] <- c("expression","sample","medium")

ggplot(pAkt1, aes(x=sample, y=expression,color=medium)) + 
  geom_boxplot() +
  stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
  theme_Publication() + 
  scale_color_manual(values=c("ES-" ="#F7C1BB", "ES+"="#F15152"))+
  labs(x="",y="z-scaled expression",title="STAT3")


