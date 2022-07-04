#GENIE3 and libraries Initialization
library(GENIE3)

#Following libraries are used for multi-core processing STILL DOESN'T WORK D:
library(doParallel)
library(foreach)
library(plyr)
library(doRNG)

set.seed(123)   #For reproducibility purposes

#Read File
path = "/media/sf_Projekt_BIONETS/federated-inference-of-grns/data/TCGA-COAD.htseq_fpkm.tsv"
OGFile = read.table(path,fileEncoding="latin1",sep="\t")


GeneT = t(t(OGFile)) #transpose matrix



#Make the format readable by GENIE 3
Dimensions = dim(GeneT)
rownum = Dimensions[1]
colnum = Dimensions[2]

rownames(GeneT) = GeneT[,1]
colnames(GeneT) = GeneT[1,]

#Eliminate "Ensemble_Id" from the matrix
GeneT[1,] = GeneT[2,]
GeneT[,1] = GeneT[,2]
GeneMat = GeneT[2:rownum,2:colnum] #Final Matrix

print('Running Genie3')

#Implementation of GENIE3
weightMat = GENIE3(GeneMat)

#matrix export 
exportpath = "/media/sf_Projekt_BIONETS/federated-inference-of-grns/data/WeightMatrix.csv"
write.table(weightMat,exportpath) 




