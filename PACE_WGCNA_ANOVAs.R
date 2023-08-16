library(afex)
library(impute)
library(dplyr)


####Set up Limma objects for DE and remove samples with missing clinical data 
R<-read.csv("Wake_Lipidyzer.csv",check.names=FALSE)
index<-apply(R[,1:10],1,is.na)
index<-apply(index,2,sum)
R<-R[index<1,]
edata<-log2(t(R[,-c(1:25)]))

####Impute missing data with KNN imputation
edata<-impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data

####Conduct SVA over samples 
pheno_table<-R[,c(2,3,14,18)]
colnames(pheno_table)<-c("StudyID","Treatment","Study","Main")
Treat<-factor(paste(pheno_table$Treatment,pheno_table$Main,sep="."))
Treat<-factor(paste(pheno_table$Study,Treat,sep="."))
Treat<-gsub("-","",Treat)
Treat<-gsub(" ","",Treat)
total<-cbind(Treat,pheno_table)

###Check for the presence of metabolites with excessively missing/low variance abundances
datExpr<-as.data.frame(t(edata))
gsg<-goodSamplesGenes(datExpr, verbose = 3)

### Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=30, by=2))
###Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


pdf(file="choice.pdf")
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

###Define soft thresholding power to be used
softPower <- 10
adjacency <- adjacency(datExpr, power = softPower)

###Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency);
dissTOM <- 1-TOM

### Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

###Define module size
minModuleSize <- 50

### Identify modules using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)

### Convert numeric labels into colors
set.seed(122)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

###Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
colors<-MEList$validColors

### Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs);

### Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average");



use<-cbind(MEs,total)
mod<-aov_car(MEblack~Treatment*Study*Main+Error(StudyID/Main),data=use)



datKME<-signedKME(datExpr, MEs, outputColumnName="MM.")

