library(GENIE3)

#Following libraries are used for multi-core processing STILL DOESN'T WORK D:
library(doParallel)
library(foreach)

library(doRNG)


export_path = "/home/hpc/iwbn/iwbn001h/GenesFRFFinal/data/Weight_Matrix.csv"

weightMat <- read.table(export_path, sep = ',')
linkList <- getLinkList(weightMat,reportMax=50000,threshold=0.1)


export_path = "/home/hpc/iwbn/iwbn001h/GenesFRFFinal/data/LinkList500.csv"

write.table(linkList, export_path, sep = ',', row.names = FALSE, col.names = FALSE)