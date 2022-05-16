dia_data <- read.table("DIA_results.csv",sep=",",header=TRUE)
dia_data <- dia_data[,-19]

convert_number <- function(x){
  x <- as.character(x)
  x <- gsub(pattern = ",", replacement = ".",x = x, fixed = TRUE)
  x <- as.numeric(x)
  return(x)
}

dia_data_number <- apply(dia_data[,c(4:18)],2,convert_number)
dia <- data.frame(dia_data$Gene, cbind(apply(dia_data_number[,1:4],1,median),
                        apply(dia_data_number[,5:8],1,median),
                        apply(dia_data_number[,9:12],1,median),
                        apply(dia_data_number[,13:15],1,median)))
colnames(dia) <- c("Gene","KO1G4", "KO3G4","LTED","MCF7")

res_KO1_LTED <- t.test(dia$KO1G4, dia$LTED, paired = FALSE, alternative = "two.sided")
