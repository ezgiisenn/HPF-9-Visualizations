library(ggplot2)
library(tidyr)

#### Read the tsv files given in html results of GSEA Analysis performed with Birgitta ###
down <- read.table("gsea_downregulated.tsv",sep="\t",header=TRUE)
up <- read.table("gsea_upregulated.tsv",sep="\t",header=TRUE)

### Filter the data based on significance ####
down <- down[,c(1,6,8)]
down <- down %>% filter(FDR.q.val < 0.25) #this FDR value can be changed based on your needs

up <- up[,c(1,6,8)]
up <- up %>% filter(FDR.q.val < 0.25)
enriched <- rbind(up,down)

### Clear the Names of Gene Sets in Hallmark###
enriched <- enriched %>% mutate(NAME = gsub("HALLMARK_","",NAME),
                    NAME=  gsub("_"," ",NAME))


pdf("enrichment_plot_from_GSEA.pdf",paper="USr")
par(mar=(c(10,5,5,5)))
ggplot(enriched,aes(x=reorder(NAME,NES),
                     y = NES,fill= FDR.q.val)) +
  geom_bar(stat="identity",width = 0.5) +
  coord_flip() +
  theme_bw() +
  labs(y = "Enrichment Scores (ES)",x="",fill = "FDR") +
  ggtitle("Downregulated Pathways        Upregulated Pathways")+
  theme(plot.title = element_text(hjust = 1,size = 15, face = "bold"),legend.position = "bottom",
        legend.direction = "horizontal",aspect.ratio = 1/2)
dev.off()
