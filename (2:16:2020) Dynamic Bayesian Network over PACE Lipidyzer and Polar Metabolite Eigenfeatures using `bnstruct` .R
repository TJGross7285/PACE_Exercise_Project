#################################(2/16/2020) Dynamic Bayesian Network over PACE Lipidyzer and Polar Metabolite Eigenfeatures using `bnstruct` 
library(dplyr) 
library(impute)
library(bnstruct)
library(WGCNA)
library(bnlearn)
##########LIPIDYZER DATA SUSPECT 5/2020
####Set up Lipidyzer data 
####setwd("~/Desktop/Results/Lipidyzer")
####index<-apply(R[,1:10],1,is.na)
####index<-apply(index,2,sum)
####R<-R[index<1,]
####index<- R$StudyID %in% c(1060202,1060127,1181043,1050103,1182023)
####R<-R[index==FALSE,]
####edata<-log2(t(R[,-c(1:25)]))
####datExpr<-t(impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data)

#### Generate final phenotypic table
E<-read.csv("Wake_Pheno.csv",check.names=FALSE)[,-c(1)]
T<-read.csv("ADD_Table.csv",check.names=FALSE)
Z<-read.csv("Baker_PACE_data.csv",check.names=FALSE,na.strings=c(".","0"))[,1:10]
Z$Sex<-as.factor(Z$Sex)
Z$APOE<-as.factor(Z$APOE)
join<-dplyr::inner_join(E,T,by="Key")
join<-dplyr::inner_join(Z,join,by="StudyID")
levels(join$APOE)<-list("HasE4"=c("2_4","3_4","4_4"),"NoE4"=c("2_3","3_3"))
levels(join$Sex)<-list("F"=c("f","F"),"M"=c("M"))
write.csv(join,file="Complete_Wake_Pheno.csv")

####Read in LC-MS data
U<-read.csv("mapstone_plasma_metabolomics_NEG.csv",check.names=FALSE,na.strings="0")
I<-read.csv("mapstone_plasma_metabolomics_POS.csv",check.names=FALSE,na.strings="0")
labelsU<-paste(paste(U[,1],U[,2],sep="_"),rep("NEG",dim(U)[1]),sep="/")
labelsI<-paste(paste(I[,1],I[,2],sep="_"),rep("POS",dim(I)[1]),sep="/")

####Format negative mode data
abunds_neg<-cbind(as.data.frame(colnames(U)[-c(1:2)]),as.data.frame(t(U[,-c(1:2)])))
colnames(abunds_neg)[1]<-"Untargeted_Sample_ID"
abunds_neg$Untargeted_Sample_ID<-sub("01$","",abunds_neg$Untargeted_Sample_ID)
join_neg<-dplyr::inner_join(join,abunds_neg,by="Untargeted_Sample_ID")
colnames(join_neg)[27:dim(join_neg)[2]]<-labelsU

####Format postive mode data
abunds_pos<-cbind(as.data.frame(colnames(I)[-c(1:2)]),as.data.frame(t(I[,-c(1:2)])))
colnames(abunds_pos)[1]<-"Untargeted_Sample_ID"
abunds_pos$Untargeted_Sample_ID<-sub("01$","",abunds_pos$Untargeted_Sample_ID)
join_pos<-dplyr::inner_join(join,abunds_pos,by="Untargeted_Sample_ID")
colnames(join_pos)[27:dim(join_pos)[2]]<-labelsI

####Integrate negative mode and positive mode data   *************FIXXXXXXXXXXXXXX!!!!!!!!!!!!!!!!!!*******************
total_frame<-cbind(join_neg,join_pos[,-c(1:26)])
index<-apply(total_frame[,3:10],1,is.na)
index<-apply(index,2,sum)
total_frame<-total_frame[index<1,]


####Generate feature metadata for negative mode 
featureData_neg<-cbind(as.data.frame(U[,1:2]),as.data.frame(rep("ESI-",dim(t(U[,1:2]))[2])))
colnames(featureData_neg)<-c("MZ","RT","Mode")

####Generate feature metadata for positive mode 
featureData_pos<-cbind(as.data.frame(I[,1:2]),as.data.frame(rep("ESI+",dim(t(I[,1:2]))[2])))
colnames(featureData_pos)<-c("MZ","RT","Mode")

####Bind feature metadata
featureData<-rbind(featureData_neg,featureData_pos)
featureData<-cbind(colnames(total_frame)[-c(1:26)],featureData)
colnames(featureData)[1]<-"FeatureID"

#####Set up metabolite abundances for DE analysis
edata<-log2(t(total_frame[,-c(1:26)]))
datExpr<-t(impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data)




















####Run WGCNA functions to cluster lipids and plot resulting modules  
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- 12
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

####Set up PACE metadata ****USE LIPIDYZER-DERIVED CODE FOR THIS!!!!!****
setwd("~/Desktop/Results/Lipidyzer")
R<-read.csv("Wake_Lipidyzer.csv",check.names=FALSE)
index<-apply(R[,1:10],1,is.na)
index<-apply(index,2,sum)
meta<-R[index<1,]


####Create combined data frame of PACE metadata and all eigenfeatures; discretize numeric features into tertiles  
final<-cbind(as.data.frame(meta),as.data.frame(All_eigenfeatures))
final_numeric<-final[,c(3:5,14,6:11,883:893)]
final_numeric$Bsln_Cholesterol<-as.numeric(final_numeric$Bsln_Cholesterol)
final_numeric$Bsln_LDL<-as.numeric(final_numeric$Bsln_LDL)
final_numeric$Bsln_HDL<-as.numeric(final_numeric$Bsln_HDL)
final_numeric$Bsln_Trigly<-as.numeric(final_numeric$Bsln_Trigly)
final_numeric$Bsln_Fasting_Glu<-as.numeric(final_numeric$Bsln_Fasting_Glu)
D<-discretize(final_numeric)

####Clean/finalize clinical metadata/discretized eigenfeatures for network modeling 
good_meta<-meta[,c(13,18)]
F<-cbind(as.data.frame(good_meta),as.data.frame(D))
colnames(F)[1:2]<-c("StudyID","Visit#")
A<-reshape(F,v.names=colnames(F)[-c(1)],timevar="Visit#",idvar="StudyID",direction="wide")
newfg<-A[,-c(1:2,24:34)]
numeric_discrete<- matrix(nrow=dim(newfg)[1],ncol=dim(newfg)[2]) 
for(i in 1:ncol(newfg)){
                numeric_discrete[,i] <- as.integer(newfg[,i])
}
numeric_discrete<-as.data.frame(numeric_discrete)

####Set up DBN components 
variables<-c("Treatment","Sex","APOE","Study",
             "Bsln_BMI","Bsln_Cholesterol","Bsln_LDL","Bsln_HDL","Bsln_Trigly","Bsln_Fasting_Glu",
             "MEblack.1","MEblue.1","MEbrown.1","MEgreen.1","MEgrey.1","MEmagenta.1","MEpink.1",    
             "MEpurple.1","MEred.1","MEturquoise.1","MEyellow.1",
             "MEblack.2","MEblue.2","MEbrown.2",
             "MEgreen.2","MEgrey.2","MEmagenta.2","MEpink.2",    
             "MEpurple.2","MEred.2","MEturquoise.2","MEyellow.2")
node<-c(rep(2,3),rep(3,29))
dataset <- BNDataset(numeric_discrete, discreteness=as.factor(rep("D",32)), variables=variables, node.sizes=node)
layers<-as.numeric(c(1,1,1,1,rep(2,17),rep(3,11)))

set.seed(122)
dataset <- bootstrap(dataset, num.boots = 1000)
dbn <- learn.network(dataset,layering=layers,scoring.func="BDeu",bootstrap=TRUE)

####Write learned network plot to PDF  
pdf(file="PACE_Layered_Network_5_29_2020.pdf")
plot(dbn)
dev.off()




















####################################################
####################################################
####################################################

####Set up polar metabolites data 
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

####Run WGCNA functions to cluster polar metabolites and plot resulting modules 
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- 9
adjacency <- adjacency(datExpr, power = softPower,type="signed")
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM
geneTree <- hclust(as.dist(dissTOM), method = "average");
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
      labels = FALSE, hang = 0.04)
minModuleSize <- 25
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                 deepSplit = 3, pamRespectsDendro = FALSE,
                 minClusterSize = minModuleSize);

####Name modules and replot; write out eigenfeature matrix 
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
MEs_Polar <- MEList$eigengenes
colorsPolar<-MEList$validColors

####Generate combined data frame of Lipidyzer/polar metabolite features 
All_eigenfeatures<-cbind(as.data.frame(MEs_Lipidyzer),as.data.frame(MEs_Polar)) 
colnames(All_eigenfeatures)<-c("MEblack_L","MEblue_L","MEbrown_L","MEgreen_L","MEgrey_L","MEpink_L","MEred_L","MEturquoise_L","MEyellow_L","MEgold_P",
                               "MEgrey_P","MEtopaz_P")

####################################################
####################################################
####################################################

####Set up PACE metadata ****USE LIPIDYZER-DERIVED CODE FOR THIS!!!!!****
setwd("~/Desktop/Results/Lipidyzer")
R<-read.csv("Wake_Lipidyzer.csv",check.names=FALSE)
index<-apply(R[,1:10],1,is.na)
index<-apply(index,2,sum)
R<-R[index<1,]
index<- R$StudyID %in% c(1060202,1060127,1181043,1050103,1182023)  
meta<-R[index==FALSE,c(2:11,14,18)]

####Create combined data frame of PACE metadata and all eigenfeatures; discretize numeric features into tertiles  
final<-cbind(as.data.frame(meta),as.data.frame(All_eigenfeatures))
final_numeric<-final[,c(5:10,13:24)]
final_numeric$Bsln_Cholesterol<-as.numeric(final_numeric$Bsln_Cholesterol)
final_numeric$Bsln_LDL<-as.numeric(final_numeric$Bsln_LDL)
final_numeric$Bsln_HDL<-as.numeric(final_numeric$Bsln_HDL)
final_numeric$Bsln_Trigly<-as.numeric(final_numeric$Bsln_Trigly)
final_numeric$Bsln_Fasting_Glu<-as.numeric(final_numeric$Bsln_Fasting_Glu)
D<-discretize(final_numeric)

####Clean/finalize clinical metadata/discretized eigenfeatures for network modeling 
good_meta<-meta[,-c(5:10)]
F<-cbind(as.data.frame(good_meta),as.data.frame(D))
A<-reshape(F,v.names=colnames(F)[-c(1,6)],timevar="Visit #",idvar="StudyID",direction="wide")
fg<-A[,-c(1)]
newfg<-fg[,-c(23:32)]
numeric_discrete<- matrix(nrow=dim(newfg)[1],ncol=dim(newfg)[2]) 
for(i in 1:ncol(newfg)){
                numeric_discrete[,i] <- as.integer(newfg[,i])
}
numeric_discrete<-as.data.frame(numeric_discrete)

####Set up DBN components 
variables<-c("Treatment","Sex","APOE","Study",
	         "Bsln_BMI","Bsln_Cholesterol","Bsln_LDL","Bsln_HDL","Bsln_Trigly","Bsln_Fasting_Glu",
	         "MEblack_L.1","MEblue_L.1",
	         "MEbrown_L.1","MEgreen_L.1",
	         "MEgrey_L.1","MEpink_L.1","MEred_L.1",
	         "MEturquoise_L.1","MEyellow_L.1","MEgold_P.1",
	         "MEgrey_P.1","MEtopaz_P.1",
	         "MEblack_L.2","MEblue_L.2",
	         "MEbrown_L.2","MEgreen_L.2",
	         "MEgrey_L.2","MEpink_L.2","MEred_L.2",
	         "MEturquoise_L.2","MEyellow_L.2","MEgold_P.2",
	         "MEgrey_P.2","MEtopaz_P.2")
node<-c(rep(2,3),rep(3,31))
dataset <- BNDataset(numeric_discrete, discreteness=as.factor(rep("D",34)), variables=variables, node.sizes=node)
layers<-as.numeric(c(1,1,1,1,rep(2,18),rep(3,12)))

set.seed(122)
dataset <- bootstrap(dataset, num.boots = 1000)
dbn <- learn.network(dataset,layering=layers,scoring.func="BDeu",bootstrap=TRUE)

####Write learned network plot to PDF  
pdf(file="PACE_Layered_Network_2_21_2020.pdf")
plot(dbn)
dev.off()


