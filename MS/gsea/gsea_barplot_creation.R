library(ggplot2)
library(tidyr)

#### Read the tsv files given in html results of GSEA Analysis performed with Birgitta ###
down <- read.table("LTED_KO1_downregulated.tsv",sep="\t",header=TRUE)
up<- read.table("LTED_KO1_gsea.tsv",sep="\t",header=TRUE)
  
### Filter the data based on significance ####
down <- down[,c(1,6,8)]
down <- down %>% filter(FDR.q.val < 0.25) #this FDR value can be changed based on your needs

up <- up[,c(1,6,8)]
up <- up %>% filter(FDR.q.val < 0.15)
enriched <- rbind(up,down)


### Clear the Names of Gene Sets in Hallmark###
enriched <- enriched %>% mutate(NAME = gsub("HALLMARK_","",NAME),
                    NAME=  gsub("_"," ",NAME))

enriched <- enriched %>% mutate(state= "down-regulated")
enriched[1:13,4] <- "up-regulated"


theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(family="Times New Roman"),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = 14),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.line = element_line(colour="#f0f0f0"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(3,3,3,3),"cm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}



ggplot(enriched,aes(x=reorder(NAME,NES),
                     y = NES,fill= state)) +
  ylim(c(-2.7,2.5)) +
  geom_bar(stat="identity",width = 0.82) +
  geom_text(data=enriched,aes(label=round(FDR.q.val,digits=4)), size = 4,
            position = position_dodge(width = 1),hjust=1.15,
            inherit.aes = TRUE,color="white") +
  theme_Publication() +
  coord_flip() +
  labs(y = "Enrichment Scores",x="",fill = "",title="LTED GLYATL1 KO1G4 vs LTED ") +
  scale_fill_manual(values=c("steelblue","tomato")) 
 

