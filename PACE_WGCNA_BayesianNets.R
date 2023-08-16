library(dplyr)
library(impute)

A<-read.csv("Complete_Wake_Pheno.csv")[,-c(1)]
B<-read.csv("20200624_pace_lipidomics data.NIST.csv",check.names=FALSE)[,-c(16)]
B<-B%>%mutate(CV = gsub("%?","",B$`Coefficient of variation`))#######%>%filter(CV<=40)

E<-read.csv("20200624_pace_lipidomics data.Abunds.csv",check.names=FALSE)
Eprime<-E%>%t()
colnames(Eprime)<-Eprime[1,]
SampleID<-colnames(E)
bound<-cbind(as.data.frame(SampleID),as.data.frame(Eprime))
abunds<-bound[2:dim(bound)[1],-c(1)]
rownames(abunds)<-NULL

F<-read.csv("20200624_pace_lipidomics data.PoolQC.csv",check.names=FALSE)[,-c(25)]
F<-F%>%mutate(CV = gsub("%?","",F$`Coefficient of variation`))########%>%filter(CV<=40)
joint<-as.character(E$compound_name) %in% intersect(as.character(B$compound_name),as.character(F$compound_name))

abunds<-abunds[,joint==TRUE]
SampleIDprime<-colnames(E)[-c(1)]

Lipidizer_Use<-cbind(as.data.frame(SampleIDprime),as.data.frame(abunds))
colnames(Lipidizer_Use)[1]<-"Key"
Lipidizer_Use<-Lipidizer_Use%>%filter(Key!="Pooled QC" & Key!="NIST Plasma")

lipid_bound<-dplyr::inner_join(A,Lipidizer_Use,by="Key")
index<-duplicated(lipid_bound$`Key`)
lipid_bound<-lipid_bound[index==FALSE,]
index<-apply(apply(lipid_bound[,1:10],2,is.na),1,sum)
lipid_bound<-lipid_bound[index<1,]


########################################
########################################
########################################
########################################


G<-read.csv("20200625_pace_metabolomics.NIST.csv",check.names=FALSE)[,-c(3,18)]
G<-G%>%mutate(CV = gsub("%?","",G$`Coefficient of variation`))%>%filter(CV<=40)

C<-read.csv("20200625_pace_metabolomics.Abunds.csv",check.names=FALSE)
Cprime<-C%>%t()
colnames(Cprime)<-Cprime[1,]
SampleID<-colnames(C)
bound<-cbind(as.data.frame(SampleID),as.data.frame(Cprime))
abunds<-bound[4:dim(bound)[1],-c(1)]
rownames(abunds)<-NULL

D<-read.csv("20200625_pace_metabolomics.PoolQC.csv",check.names=FALSE)[,-c(3,18)]
D<-D%>%mutate(CV = gsub("%?","",D$`Coefficient of variation`))%>%filter(CV<=40)
joint<-as.character(C$`Unique names`) %in% intersect(as.character(G$`Unique names`),as.character(D$`Unique names`))

abunds<-abunds[,joint==TRUE]
SampleIDprime<-colnames(C)[-c(1:3)]

Metabolites_Use<-cbind(as.data.frame(SampleIDprime),as.data.frame(abunds))
colnames(Metabolites_Use)[1]<-"Key"
Metabolites_Use<-Metabolites_Use%>%filter(Key!="Pooled QC" & Key!="NIST Plasma")

metab_bound<-dplyr::inner_join(A,Metabolites_Use,by="Key")
index<-duplicated(metab_bound$`Key`)
metab_bound<-metab_bound[index==FALSE,]
index<-apply(apply(metab_bound[,1:10],2,is.na),1,sum)
metab_bound<-metab_bound[index<1,]


all.equal(metab_bound$StudyID,lipid_bound$StudyID)
write.csv(lipid_bound,file="Clean_PACE_Lipid_9_23.csv")
write.csv(metab_bound,file="Clean_PACE_Metab_9_23.csv")



A<-read.csv("Complete_Wake_Pheno.csv")[,-c(1)]
W<-read.csv("PACE_Proteomics_7_2020.csv",check.names=FALSE)
########U<-scale(log2(W[,-c(1:4)]))
U<-cbind(W[,1:4],as.data.frame(W[,-c(1:4)]))
join.fac<-as.factor(paste(W$StudyID,W$`Visit`,sep="."))
U<-cbind(join.fac,U)

join.fac<-as.factor(paste(A$StudyID,A$`Visit..`,sep="."))
A<-cbind(join.fac,A)

prot_bound<-dplyr::inner_join(A,U,by="join.fac")
index<-duplicated(prot_bound$`Key`)
prot_bound<-prot_bound[index==FALSE,]
index<-apply(apply(prot_bound[,5:11],2,is.na),1,sum)
prot_bound<-prot_bound[index<1,]
write.csv(prot_bound,file="Clean_PACE_Prot_9_23.csv")

three.way<-as.factor(paste(paste(prot_bound$Study,prot_bound$Treatment,sep="."),prot_bound$`Visit..`,sep="."))
abunds<-prot_bound[,-c(1:31)]

pdf(file="PACE_Proteomics_VIO_7_25.pdf")
for (i in 1:dim(abunds)[2]){
	vioplot(abunds[,i]~three.way, main=colnames(abunds)[i])
}
dev.off()

mat <- matrix(, nrow = dim(abunds)[2], ncol = 1)
for(i in 1:dim(abunds)[2]){
  mat[i,]<- bd.test(abunds[,i],y=as.numeric(three.way))$p.value
}
mat<-cbind(as.data.frame(colnames(abunds)),as.data.frame(mat))




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
final<-cbind(as.data.frame(pheno),as.data.frame(All_eigenfeatures))
colnames(final)[18]<-"Main"



























A<-read.csv("Complete_Wake_Pheno.csv")[,-c(1)]
E<-read.csv("20200624_pace_lipidomics data.Abunds.csv",check.names=FALSE,header=FALSE,na.strings=c("","NA","."))
abunds<-as.data.frame(t(E[-c(1),-c(1)]))
colnames(abunds)<-E[-c(1),1]
Key<-as.vector(unlist(E[1,]))[-c(1)]

bound<-cbind(as.data.frame(Key),as.data.frame(abunds))%>%filter(Key!="NIST Plasma" & Key!="Pooled QC")
metab_bound<-dplyr::inner_join(A,bound,by="Key")
index<-duplicated(metab_bound$`Tube..`)

metab_bound<-metab_bound[index==FALSE,]
write.csv(metab_bound,file="No_CV_LIP.csv")























A<-read.csv("Clean_PACE_Lipid_7_14.csv",check.names=FALSE)
edata<-log2(t(A[,-c(1:27)]))
lipid_datExpr<-scale(t(impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data))



B<-read.csv("Clean_PACE_Metab_7_14.csv",check.names=FALSE)
edata<-log2(t(B[,-c(1:27)]))
metab_datExpr<-scale(t(impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data))

datExpr<-cbind(lipid_datExpr,metab_datExpr)

####Run WGCNA functions to cluster lipids and plot resulting modules  
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- 3
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
final<-cbind(as.data.frame(A[,2:27]),as.data.frame(All_eigenfeatures))
colnames(final)[17]<-"Main"










##############################
##############################
##############################
##############################
##############################


A<-read.csv("Complete_Wake_Pheno.csv")[,-c(1)]
E<-read.csv("20200624_pace_lipidomics data.Abunds.csv",check.names=FALSE,header=FALSE,na.strings=c("","NA","."))
abunds<-as.data.frame(t(E[-c(1),-c(1)]))
colnames(abunds)<-E[-c(1),1]
Key<-as.vector(unlist(E[1,]))[-c(1)]

bound<-cbind(as.data.frame(Key),as.data.frame(abunds))%>%filter(Key!="NIST Plasma" & Key!="Pooled QC")
metab_bound<-dplyr::inner_join(A,bound,by="Key")
index<-duplicated(metab_bound$`Tube..`)

metab_bound<-metab_bound[index==FALSE,]
write.csv(metab_bound,file="No_CV_LIP.csv")








library(dplyr)
library(WGCNA)
library(impute)
library(bnlearn)
library(bnstruct)

R<-read.csv("No_CV_LIP.csv",check.names=FALSE)
index<-apply(apply(R[,1:27],1,is.na),2,sum)
R<-R[index<1,]
edata<-log2(t(R[,-c(1:27)]))
datExpr<-t(impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data)



####Run WGCNA functions to cluster lipids and plot resulting modules  
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- 10
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


final<-cbind(as.data.frame(R[,1:27]),as.data.frame(All_eigenfeatures))
colnames(final)[18]<-"Main"





####Create combined data frame of PACE metadata and all eigenfeatures; discretize numeric features into tertiles  
final_numeric<-final[,c(3:11,14,28:29,31:32)]
final_numeric$Bsln_BMI<-as.numeric(final_numeric$Bsln_BMI)
final_numeric$Bsln_Cholesterol<-as.numeric(final_numeric$Bsln_Cholesterol)
final_numeric$Bsln_LDL<-as.numeric(final_numeric$Bsln_LDL)
final_numeric$Bsln_HDL<-as.numeric(final_numeric$Bsln_HDL)
final_numeric$Bsln_Trigly<-as.numeric(final_numeric$Bsln_Trigly)
final_numeric$Bsln_Fasting_Glu<-as.numeric(final_numeric$Bsln_Fasting_Glu)
S<-apply(final_numeric[,-c(1:3,10)],2,as.numeric)%>%as.data.frame()
D<-cbind(as.data.frame(final$Main),as.data.frame(final$StudyID),as.data.frame(final$Study),as.data.frame(final_numeric[,c(1:3,10)]),as.data.frame(discretize(S)))
colnames(D)[1:3]<-c("Visit #","StudyID","Study")
A<-reshape(D,v.names=colnames(D)[-c(1:2)],timevar="Visit #",idvar="StudyID",direction="wide")

newfg<-A[,-c(1:2,13:15,17:27)]
numeric_discrete<- matrix(nrow=dim(newfg)[1],ncol=dim(newfg)[2]) 
for(i in 1:ncol(newfg)){
                numeric_discrete[,i] <- as.integer(newfg[,i])
}
numeric_discrete<-as.data.frame(numeric_discrete)

####Set up DBN components 
variables<-c("Treatment","Sex","APOE","Study",
	         "Bsln_BMI","Bsln_Cholesterol","Bsln_LDL","Bsln_HDL","Bsln_Trigly","Bsln_Fasting_Glu","MEyellow.Pre",
             "MEblue.Post","MEbrown.Post","MEturquoise.Post","MEyellow.Post")
node<-c(rep(2,3),3,rep(3,11))
dataset <- BNDataset(numeric_discrete, discreteness=as.factor(rep("D",15)), variables=variables, node.sizes=node)
layers<-as.numeric(c(1,1,1,1,rep(2,7),rep(3,4)))

set.seed(122)
dataset <- bootstrap(dataset, num.boots = 1000)
dbn <- learn.network(dataset,layering=layers,scoring.func="BDeu",bootstrap=TRUE)

####Write learned network plot to PDF  
pdf(file="PACE_JHKH_7_10.pdf")
plot(dbn)
dev.off()
