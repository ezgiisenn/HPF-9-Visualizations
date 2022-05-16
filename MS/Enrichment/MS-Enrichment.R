library(readxl)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(pheatmap)
library(scales)
library(fgsea)
library(biomaRt)
library(msigdbr)
library(clusterProfiler)
library(tibble)



KO1_LTED <- read.delim("significant_KO_1G4_vs_LTED.txt",sep="\t",header=TRUE)
KO2_LTED <- read.delim("significant_KO_3G4_vs_LTED.txt",sep="\t",header=TRUE)
LTED_WT <- read.delim("significant_LTED_vs_WT.txt",sep="\t",header=TRUE)

convert_number <- function(x){
  x <- as.character(x)
  x <- gsub(pattern = ",", replacement = ".",x = x, fixed = TRUE)
  x <- as.numeric(x)
  return(x)
}

data_adjust <- function(data){
  data <- data.frame(cbind(data[,24],data[,21],data[,20]))
  colnames(data) <- c("gene_name","log2fc","p_val")
  data[,2] <- convert_number(data[,2])
  data[,3] <- convert_number(data[,3])
  data <- data %>% filter(p_val < 0.05) %>%
              arrange(desc(log2fc))
  return(data)
}

write.table(data_adjust(KO1_LTED)[,1:2],"KO1_LTED.txt",sep="\t",row.names=FALSE,quote=FALSE)
write.table(data_adjust(KO2_LTED)[,1:2],"KO2_LTED.txt",sep="\t",row.names=FALSE,quote=FALSE)
write.table(data_adjust(LTED_WT)[,1:2],"LTED_WT.txt",sep="\t",row.names=FALSE,quote=FALSE)


"""
####### Knock out vs LTED (LogFold indicates Knockout's Expression) ########

KO1_LTED_data <- data.frame(cbind(KO1_LTED$T..Gene,KO1_LTED$N...Log.Student.s.T.test.p.value.KO_1G4_LTED,
                                  KO1_LTED$N..Student.s.T.test.Difference.KO_1G4_LTED))
colnames(KO1_LTED_data) <- c("gene_name","p_val","log2fc")


KO1_LTED_data[,2] <- convert_number(KO1_LTED_data[,2])
KO1_LTED_data[,3] <- convert_number(KO1_LTED_data[,3])

symbols <- mapIds(org.Hs.eg.db, keys = KO1_LTED_data$gene_name, keytype = "SYMBOL", column="ENSEMBL")
KO1_LTED_data <- cbind(KO1_LTED_data,data.frame(symbols))
KO1_LTED_data <- KO1_LTED_data[-which(KO1_LTED_data$gene_name =="TWISTNB"),]

KO1_LTED_data <- KO1_LTED_data %>% arrange(desc(log2fc))

H <-  msigdbr(species="Homo sapiens",category="H")

H.ensembl.ls <- H %>%
  select(gs_name,ensembl_gene) %>%
  group_by(gs_name) %>%
  summarise(zll.genes = list(unique(ensembl_gene))) %>%
  deframe()

GSEA.vec <- KO1_LTED_data$log2fc
names(GSEA.vec) <- KO1_LTED_data$symbols

#min(GSEA.vec)
#max(GSEA.vec)

gsea.H <- fgseaSimple(pathways = H.ensembl.ls,
                      stats = GSEA.vec,
                      scoreType = "std",
                      nperm = 1000)

data_gsea <- gsea.H %>%
  filter(padj <= 0.8) %>%
  mutate(pathway = gsub("HALLMARK_","",pathway),
         pathway = gsub("_"," ",pathway))

#pdf("enrichment_plot_R.pdf",paper="USr")
par(mar=(c(10,10,5,5)))
ggplot(data_gsea,aes(x=reorder(pathway,ES),
                     y = ES,fill= padj)) +
  geom_bar(stat="identity",width = 0.5) +
  coord_flip() +
  theme_bw() +
  labs(y = "Enrichment Scores (ES)",x="",fill = "FDR") +
  ggtitle("Pathways in KO2 vs LTED ")+
  theme(plot.title = element_text(hjust = 0,size = 15, face = "bold"),legend.position = "bottom",
        legend.direction = "horizontal",aspect.ratio = 1/2)
dev.off()




########## 2. LTED vs WT ###################

LTED_WT_data <- data.frame(cbind(LTED_WT$T..Gene,LTED_WT$N...Log.Student.s.T.test.p.value.LTED_WT,
                                 LTED_WT$N..Student.s.T.test.Difference.LTED_WT))


colnames(LTED_WT_data) <- c("gene_name","p_val","log2fc")

LTED_WT_data[,2] <- convert_number(LTED_WT_data[,2])
LTED_WT_data[,3] <- convert_number(LTED_WT_data[,3])

symbols <- mapIds(org.Hs.eg.db, keys = LTED_WT_data$gene_name, keytype = "SYMBOL", column="ENSEMBL")
LTED_WT_data <- cbind(LTED_WT_data,data.frame(symbols))

LTED_WT_data <- LTED_WT_data[-c(8,41,154),]
LTED_WT_data <- LTED_WT_data %>% arrange(desc(log2fc))

write.table(LTED_WT_data$symbols,"LTED_WT.txt",row.names=FALSE,quote=FALSE)
write.table(KO1_LTED_data$symbols,"KO1_LTED.txt",row.names=FALSE,quote=FALSE)

H <-  msigdbr(species="Homo sapiens",category="H")

H.ensembl.ls <- H %>%
  select(gs_name,ensembl_gene) %>%
  group_by(gs_name) %>%
  summarise(zll.genes = list(unique(ensembl_gene))) %>%
  deframe()

indx <- which(LTED_WT_data$symbols == unique(LTED_WT_data$symbols))
LTED_WT_data <- LTED_WT_data[indx,]
GSEA.vec <- LTED_WT_data$log2fc
names(GSEA.vec) <- LTED_WT_data$symbols

####### 3  KO2 vs LTED###########

KO2_LTED_data <- data.frame(cbind(KO2_LTED$T..Gene,KO2_LTED$N...Log.Student.s.T.test.p.value.KO_3G4_LTED,
                       KO2_LTED$N..Student.s.T.test.Difference.KO_3G4_LTED))
colnames(KO2_LTED_data) <- c("gene_name","p_val","log2fc")


KO2_LTED_data[,2] <- convert_number(KO2_LTED_data[,2])
KO2_LTED_data[,3] <- convert_number(KO2_LTED_data[,3])

symbols <- mapIds(org.Hs.eg.db, keys = KO2_LTED_data$gene_name, keytype = "SYMBOL", column="ENSEMBL")
KO2_LTED_data <- cbind(KO2_LTED_data,data.frame(symbols))

KO2_LTED_data <- KO2_LTED_data[-c(18,43),]

GSEA.vec <- KO2_LTED_data$log2fc
names(GSEA.vec) <- KO2_LTED_data$symbols

indx <- which(KO2_LTED_data$symbols == unique(KO2_LTED_data$symbols))
KO2_LTED_data <- KO2_LTED_data[indx,]



"""






