#######################Set Up Lipidyzer WGCNA 
setwd("~/Desktop/Results/Lipidyzer")
R<-read.csv("Wake_Lipidyzer.csv",check.names=FALSE)
index<-apply(R[,1:10],1,is.na)
index<-apply(index,2,sum)
R<-R[index<1,]
index<- R$StudyID %in% c(1060202,1060127,1181043,1050103,1182023)
R<-R[index==TRUE,]
edata<-log2(t(R[,-c(1:25)]))
datExpr<-t(impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data)



setwd("~/Desktop/Results")
R<-read.csv("Energy_PACE_Fixed.csv",check.names=FALSE)
index<-apply(R[,1:11],1,is.na)
index<-apply(index,2,sum)
R<-R[index<1,]
thresh_index<-apply(apply(R,2,is.na),2,sum)/dim(R)[1]
R<-R[,thresh_index<.333]
index<- R$StudyID %in% c(1060202,1060127,1181043,1050103,1182023)
R<-R[index==FALSE,]
edata<-log2(t(R[,-c(1:29)]))
datExpr<-t(impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data)


#######################Set Up Lipidyzer WGCNA 
setwd("~/Desktop/Results/Lipidyzer")
R<-read.csv("Wake_Lipidyzer.csv",check.names=FALSE)
index<-apply(R[,1:10],1,is.na)
index<-apply(index,2,sum)
R<-R[index<1,]
index<- R$StudyID %in% c(1060202,1060127,1181043,1050103,1182023)  
meta<-R[index==TRUE,]