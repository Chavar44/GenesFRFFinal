#GENIE3 and libraries Initialization
library(GENIE3)

#Following libraries are used for multi-core processing STILL DOESN'T WORK D:
library(doParallel)
library(foreach)
library(plyr)
library(doRNG)

set.seed(123)   #For reproducibility purposes

data_path <- "/media/sf_Projekt_BIONETS/federated-inference-of-grns/data"


#Read File
path <- paste0(data_path, "/TCGA-COAD.htseq_fpkm.tsv")
data <- read.table(path, fileEncoding = "latin1", sep = "\t")

#transpose matrix
data <- t(t(data))


#Make the format readable by GENIE 3
Dimensions <- dim(data)
rownum <- Dimensions[1]
colnum <- Dimensions[2]

data[1:rownum, 1] = substr(data[1:rownum], 1, 15) #Remove decimal points from Genes Ensembl id

rownames(data) <- data[, 1]
colnames(data) <- data[1,]

#Eliminate "Ensemble_Id" from the matrix
data[1,] <- data[2,]
data[, 1] <- data[, 2]
data <- data[2:rownum, 2:colnum] #Final Matrix

#Read Regulators
path <- paste0(data_path, "/Regulators.txt")
Regulators <- read.table(path, fileEncoding = "latin1", sep = "\n")
Regulators <- Regulators[, 1]

#Implementation of GENIE3
weightMat <- GENIE3(data, regulators = Regulators, verbose=TRUE)

#matrix export
export_path <- paste0(data_path, "/Weight_Matrix.csv")
write.table(weightMat, exportpath, sep = ',', row.names = FALSE, col.names = FALSE)
cat(0)




