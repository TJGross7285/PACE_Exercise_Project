library(impute)
library(WGCNA)
library(afex)


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

merge<-cbind(as.data.frame(pheno),as.data.frame(final[,32:51]),as.data.frame(All_eigenfeatures))
colnames(merge)[18]<-"Main"

########EXCLUDE IL-5 need to fix 
discretized<-discretize(as.data.frame(apply(merge[,c(6:11,32:40,42:43,45:53)],2,as.numeric)))%>%as.data.frame()




D<-cbind(as.data.frame(merge$Main),as.data.frame(merge$StudyID.x),as.data.frame(merge$Study),as.data.frame(merge[,3:5]),as.data.frame(discretized))
colnames(D)[1:3]<-c("Visit #","StudyID","Study")
A<-reshape(D,v.names=colnames(D)[-c(1:2)],timevar="Visit #",idvar="StudyID",direction="wide")

newfg<-A[,-c(1:5,32:41)]
numeric_discrete<- matrix(nrow=dim(newfg)[1],ncol=dim(newfg)[2]) 
for(i in 1:ncol(newfg)){
                numeric_discrete[,i] <- as.integer(newfg[,i])
}
numeric_discrete<-as.data.frame(numeric_discrete)


addnl<-cbind(as.data.frame(A[,3:5]),as.data.frame(A$`Study.1`))
colnames(addnl)[4]<-"Study.1"
addnl$Treatment.1<-as.factor(addnl$Treatment.1)
addnl$APOE.1<-as.factor(addnl$APOE.1)
addnl$Study.1<-as.factor(addnl$Study.1)
addnl$Sex.1<-as.factor(addnl$Sex.1)

addnl$Treatment.1<-as.numeric(addnl$Treatment.1)
addnl$APOE.1<-as.numeric(addnl$APOE.1)
addnl$Study.1<-as.numeric(addnl$Study.1)
addnl$Sex.1<-as.numeric(addnl$Sex.1)

final_discrete<-cbind(addnl,numeric_discrete)
####Set up DBN components 
variables<-c("Treatment","Sex","APOE","Study",
	         "Bsln_BMI","Bsln_Cholesterol","Bsln_LDL","Bsln_HDL","Bsln_Trigly","Bsln_Fasting_Glu",
	         paste(colnames(as.data.frame(final[,32:51])[-c(10,13)]),"1",sep="."),"MEblue.1","MEturquoise.1",
	         paste(colnames(as.data.frame(final[,32:51])[-c(10,13)]),"2",sep="."), "MEblue.2","MEturquoise.2")
node<-c(rep(2,3),3,rep(3,26),rep(3,20))
dataset <- BNDataset(final_discrete, discreteness=as.factor(rep("D",50)), variables=as.character(variables), node.sizes=node)
layers<-as.numeric(c(1,1,1,1,rep(2,26),rep(3,20)))

set.seed(122)
dataset <- bootstrap(dataset, num.boots = 1000)
dbn <- learn.network(dataset,layering=layers,scoring.func="BDeu",bootstrap=TRUE)














igt<-final%>%filter(Study=="PACE-IGT")
mci<-final%>%filter(Study=="PACE-MCI")
too<-final%>%filter(Study=="PACE 2")

test<-aov_car(MEblue~Treatment*Main+Error(StudyID.x/Main),data=igt)

datExpr[,c(1:20)]



