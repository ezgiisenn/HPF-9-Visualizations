library(readxl)
library(pheatmap)
library(tidyr)

data <- read_excel("cool_practical-NT_median_normdatXexpression_log2.xlsx")
data <- data.frame(data[,-1])

data <- data %>% mutate(sample_name= gsub("_"," ",sample_name))

rownames(data) <- data[,1]
data <- data.frame(data[,-1])

data <- data.frame(data[-c(1,4,7,10),])
data <- data.frame(data[,-c(11,12)])

pdf("rppa_heatmap_without_sample1",paper="USr")
par(mar = c(10, 10, 10, 10))
pheatmap(t(data))
dev.off()
