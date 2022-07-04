#GENIE3 and libraries Initialization
library(GENIE3)

#Following libraries are used for multi-core processing STILL DOESN'T WORK D:
library(doParallel)
library(foreach)
library(plyr)
library(doRNG)

set.seed(123)   #For reproducibility purposes

#Read File
path = "C:/HMDA/Proyecto Random Forest/repository/federated-inference-of-grns/genie3/TCGA-COAD.htseq_fpkm.tsv"
OGFile = read.table(path,fileEncoding="latin1",sep="\t")

 #transpose matrix
GeneT = t(OGFile)
GeneT = t(GeneT)




#Make the format readable by GENIE 3
Dimensions = dim(GeneT)
rownum = Dimensions[1]
colnum = Dimensions[2]

GeneT[1:rownum,1] = substr(GeneT[1:rownum],1,15) #Remove decimal points from Genes Ensembl id

rownames(GeneT) = GeneT[,1]
colnames(GeneT) = GeneT[1,]

#Eliminate "Ensemble_Id" from the matrix
GeneT[1,] = GeneT[2,]
GeneT[,1] = GeneT[,2]
GeneMat = GeneT[2:rownum,2:colnum] #Final Matrix

#Read Regulators
path = "C:/HMDA/Proyecto Random Forest/repository/federated-inference-of-grns/genie3/Regulators.txt"
Regulators = read.table(path,fileEncoding="latin1",sep="\n")
Regulators = Regulators[,1]

#Implementation of GENIE3
weightMat = GENIE3(GeneMat,verbose=TRUE, regulators = Regulators)

#matrix export
#exportpath = "C:/HMDA/Proyecto Random Forest/repository/federated-inference-of-grns/src/evaluation/Data.csv"
#write.csv(Data,exportpath,row.names=FALSE)
exportpath = "C:/HMDA/Proyecto Random Forest/repository/federated-inference-of-grns/src/evaluation/WeightMatrix.csv"
write.csv(weightMat,exportpath,row.names=FALSE)




