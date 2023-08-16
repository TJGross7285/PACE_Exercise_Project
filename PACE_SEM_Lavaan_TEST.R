######################Levaan 
#############recode as one hot binary categoricals 
library(lavaan)
library(dplyr)
library(WGCNA)
library(impute)
library(bnlearn)
library(bnstruct)

A<-read.csv("Clean_PACE_Lipid_8_11.csv",check.names=FALSE)
edata<-log2(t(A[,-c(1:27)]))
lipid_datExpr<-scale(t(impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data))



B<-read.csv("Clean_PACE_Metab_8_11.csv",check.names=FALSE)
edata<-log2(t(B[,-c(1:27)]))
metab_datExpr<-scale(t(impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data))

prot_datExpr<-read.csv("Clean_PACE_Prot_8_4.csv",check.names=FALSE)
colnames(prot_datExpr)[2]<-"index"
prot_datExpr$index<-as.factor(prot_datExpr$index)

index<-paste(A$StudyID,A$`Visit..`,sep=".")
datExpr<-cbind(as.data.frame(index),as.data.frame(lipid_datExpr),as.data.frame(metab_datExpr))

final<-dplyr::inner_join(prot_datExpr[,-c(1)],datExpr,by="index")

datExpr<-final[,-c(1:31)]
pheno<-final[,1:31]

####for lip and metab 
datExpr<-datExpr[,-c(1:20)]


####Run WGCNA functions to cluster lipids and plot resulting modules  
powers <- c(c(1:10), seq(from = 12, to=30, by=2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- 4
adjacency <- adjacency(datExpr, power = softPower,type="signed")
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
minModuleSize <- 50
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                 deepSplit = 2, pamRespectsDendro = FALSE,
                 minClusterSize = minModuleSize);

####Name modules and replot; write out eigenlipid matrix 
table(dynamicMods)
set.seed(122)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                     dendroLabels = FALSE, hang = 0.03,
                     addGuide = TRUE, guideHang = 0.05,
                     main = "Gene dendrogram and module colors")
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
All_eigenfeatures <- MEList$eigengenes
colorsUntargeted<-MEList$validColors

all.equal(A$SampleID,B$SampleID)
final<-cbind(as.data.frame(pheno[,c(18,2,14,3:11)]),as.data.frame(All_eigenfeatures))
colnames(final)[1:3]<-c("Visit #","StudyID","Study")

A<-reshape(final,v.names=colnames(D)[-c(1:2)],timevar="Visit #",idvar="StudyID",direction="wide")

use<-A[,c(1:13,24:25)]

model<-as.data.frame(model.matrix(~as.factor(Study.1)+as.factor(Treatment.1)+as.factor(Sex.1)
					+as.factor(APOE.1)+Bsln_BMI.1+Bsln_Cholesterol.1+Bsln_LDL.1+Bsln_HDL.1+Bsln_Trigly.1
					+Bsln_Fasting_Glu.1+MEblue.1+MEturquoise.1+MEblue.2+MEturquoise.2,data=use))
colnames(model)[1:6]<-c("Intercept","PACE-IGT","PACE-MCI","Stretch","Male","NoE4")

model[,7:16]<-scale(model[,7:16])



sem_model<-
"i =~ 1*`MEblue.1` + 1*`MEblue.2` 
s =~ 0*`MEblue.1` + 1*`MEblue.2` 

i ~ `PACE-IGT`+`PACE-MCI`+`Stretch`
s ~ `PACE-IGT`+`PACE-MCI`+`Stretch`"

fit <- growth(sem_model, data=model)

