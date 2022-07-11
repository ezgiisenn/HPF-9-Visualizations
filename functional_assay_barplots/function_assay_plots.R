library(readxl)
library(dplyr)
library(statVisual)
library(ggthemes)
library(ggpattern)
"""
functional_data <- read_excel("Total_Cell_Counts.xlsx",col_names=FALSE)
ES_neg <- functional_data[1:6,]
day2 <- ES_neg[,3]
ES_neg <- ES_neg[,-3]
ES_neg <- rbind(ES_neg,ES_neg)
ES_neg[7:12,2] <- day2
ES_neg <- cbind(ES_neg,c("1","1","1","1","1","1","2","2","2","2","2","2"))
colnames(ES_neg) <- c("sample","day_value","day_annot")
"""

ES_neg <- read_excel("pi_divided_by_total_ES-median.xlsx",col_names=FALSE)
colnames(ES_neg) <- c("sample","day_value","day_annot")

ggplot(data=ES_neg,aes(x=day_annot,y=day_value,color=sample)) +
  geom_point(shape=15,size=5) +
  geom_line(aes(group = sample))+
  labs(x = "day", y = "total cell number") +
  scale_color_manual(values = c("MCF7 ES-" = "tan1", 
                                "LTED ES-" = "seagreen",
                                "LTED KO 1G4 ES-" = "pink",
                                "LTED KO 3G4 ES-" = "purple",
                                "LTED KD si ctrl ES-" = "lightblue",
                                "LTED KD ES-" = "darkblue",
                                "MCF7 ES-" = "tan1", 
                                "LTED ES-" = "seagreen",
                                "LTED KO 1G4 ES-" = "pink",
                                "LTED KO 3G4 ES-" = "purple",
                                "LTED KD si ctrl ES-" = "lightblue",
                                "LTED KD ES-" = "darkblue")) +
  theme_Publication(14,"sans") +
  theme(
    plot.margin=unit(c(10, 10, 5, 10), units="line"),
    legend.margin=unit(0, "lines"))

ES_pos <- functional_data[7:12,]
day2 <- ES_pos[,3]
ES_pos <- ES_pos[,-3]
ES_pos <- rbind(ES_pos,ES_pos)
ES_pos[7:12,2] <- day2
ES_pos <- cbind(ES_pos,c("1","1","1","1","1","1","2","2","2","2","2","2"))
colnames(ES_pos) <- c("sample","day_value","day_annot")

ggplot(data=ES_pos,aes(x=day_annot,y=day_value,color=sample)) +
  geom_point(shape=15,size=5) +
  geom_line(aes(group = sample))+
  labs(x = "day", y = "total cell number") +
  scale_color_manual(values = c("MCF7 ES+" = "tan1", 
                                "LTED ES+" = "seagreen",
                                "LTED KO 1G4 ES+" = "pink",
                                "LTED KO 3G4 ES+" = "purple",
                                "LTED KD si ctrl ES+" = "lightblue",
                                "LTED KD ES+" = "darkblue",
                                "MCF7 ES+" = "tan1", 
                                "LTED ES+" = "seagreen",
                                "LTED KO 1G4 ES+" = "pink",
                                "LTED KO 3G4 ES+" = "purple",
                                "LTED KD si ctrl ES+" = "lightblue",
                                "LTED KD ES+" = "darkblue")) +
  theme_Publication(14,"sans") +
  theme(
    plot.margin=unit(c(10, 10, 5, 10), units="line"),
    legend.margin=unit(0, "lines"))


ES_neg <- data.frame(functional_data[1:6,])
colnames(ES_neg) <- c("sample","day1","day2")

Es_neg_Day1 <- ES_neg[,1:2]
Es_neg_Day2 <- cbind(ES_neg[,1],ES_neg[,3])

ES_neg_new <- data.frame(rbind(data.frame(cbind(ES_neg[,1],ES_neg[,2])),
                               data.frame(cbind(ES_neg[,1],ES_neg[,3]))))
ES_neg_new[1:6,3] <- "day 1"
ES_neg_new[7:12,3] <- "day 2"                         
colnames(ES_neg_new) <- c("sample","cell_no","day")
ES_neg_new$cell_no <- as.numeric(ES_neg_new$cell_no)
"""
theme_Publication <- function(base_size=14, base_family="Times New Roman") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = "darkgray"),
            axis.title.y = element_text(angle=90,vjust =1),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x= element_text(face="bold"),
            axis.text = element_text(size=12), 
            axis.line = element_line(colour="black",size=0.4),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.position = "bottom",
            legend.title = element_blank(),
            plot.margin=unit(c(1,2,2,2),"cm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}



read_data <- function(name){
  functional_data <- data.frame(read_excel(paste(name,".xlsx",sep=""),col_names=FALSE))
  colnames(functional_data) <- c("sample","cell_no","stdev","day")
  functional_data$cell_no <- as.numeric(functional_data$cell_no)
  functional_data$stdev <- as.numeric(functional_data$stdev)
  return(functional_data)
}


plot_bar <- function(data){
  p <- ggplot(data, aes(x=reorder(sample,-cell_no),y=cell_no,fill=sample,pattern=day)) +
  geom_bar_pattern(stat="identity",
                   position ="dodge",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6,) +
  geom_errorbar(aes(ymin=cell_no-stdev, ymax=cell_no+stdev), width=0.25,size=1,
                position=position_dodge(0.9)) +
  theme_Publication() +
  scale_pattern_manual(values = c("day 1" = "circle" ,"day 4" = "stripe", "day 7" = "none")) + 
  scale_fill_manual(values=c("LTED ES+" = "#F9DC5C",
                             "LTED KD ES+" = "#087E8B",
                             "LTED KD si ctrl ES+"= "#087E8B",
                             "LTED KO 1G4 ES+"= "#EF7B45",
                             "LTED KO 3G4 ES+"= "#EF7B45",
                             "MCF7 ES+"= "#9E788F")) +
  labs(x="Apoptotic Cells in ES+ DMEM",y=expression(atop("Apoptotic Cell ", 
                                    "Percentage (%)")), title= "Apoptotic Cell Percentage in ES+ DMEM") 
  #theme(axis.text.x = element_text(angle = 60,vjust = 1, hjust=1)) 
  #guides(pattern = guide_legend(override.aes = list(fill = "white")),
         #fill = guide_legend(override.aes = list(pattern = "none")))
  return(p)
}


p1 <- plot_bar(read_data("proliferation_ES_neg"))
p2 <- plot_bar(read_data("proliferation_ES_pos"))

data <- read_data("pi_divided_by_total_ES+median")
data <- data[-which(data$day == "day 0"),]
p3 <- plot_bar(data)
p4 <- plot_bar(read_data("pi_divided_by_total_ES+median"))


pi_plot_grid(pi_ES_neg,pi_ES_pos,align="h",rel_widths=c(0.1,0.1))

library(stringr)
#Across Cells Proliferation
cells <- c(6.03393167,6.35534212,15.9913051,35.6853353)
cells <- as.numeric(cells)
medium <- c("ES-","ES-","ES+","ES+")
days <- c("day 4", "day 7", "day 4", "day 7")
stdv <- c(0.43,0.46,2.07,4.66)
data <- data.frame(cbind(medium,cells,days,stdv))
data$cells = as.double(data$cells) # <-- converting 
data$stdv = as.double(data$stdv) # <-- converting 

ggplot(data,aes(x=medium,y=cells,fill=medium,pattern=days)) + 
  geom_bar_pattern(stat="identity",position=position_dodge(0.45),width=0.45,
                   color = "white", 
                   pattern_fill = "lightgrey",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6,) +
  scale_fill_manual(values=c("ES-"="#F7C1BB","ES+"="#F15152")) +
  theme_Publication() +
  geom_errorbar(aes(ymin=cells-stdv, ymax=cells+stdv), width=0.3,
                position=position_dodge(0.5)) +
  scale_pattern_manual(values = c("day 4" = "stripe","day 7" = "none")) +
  labs(x="",title="Proliferation Across all Samples",
        y=str_pad(expression("Cell Number \n (Relative to Seed Control)"),side="right",width=0.2))

#Apoptosis Across All Cells       
cells <- c(27.53,37.83,22.11,28.49)
cells <- as.numeric(cells)
medium <- c("ES-","ES-","ES+","ES+")
days <- c("day 4", "day 7", "day 4", "day 7")
stdv <- c(2.75,5.14,1.66,1.53)
data <- data.frame(cbind(medium,cells,days,stdv))
data$cells = as.double(data$cells) # <-- converting 
data$stdv = as.double(data$stdv) # <-- converting 

ggplot(data,aes(x=medium,y=cells,fill=medium,pattern=days)) + 
  geom_bar_pattern(stat="identity",position=position_dodge(0.45),width=0.45,
                   color = "white", 
                   pattern_fill = "lightgrey",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6,) +
  geom_errorbar(aes(ymin=cells-stdv, ymax=cells+stdv), width=0.3,
                position=position_dodge(0.5)) +
  scale_fill_manual(values=c("ES-"="#F7C1BB","ES+"="#F15152")) +
  theme_Publication() +
  scale_pattern_manual(values = c("day 4" = "stripe","day 7" = "none")) +
  labs(x="",title="Apoptosis Across all Samples",y=str_pad(expression("Cell Number \n (Relative to Seed Control)"),side="right",width=0.2))


