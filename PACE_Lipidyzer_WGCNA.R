#######################Set Up Lipidyzer WGCNA 
setwd("~/Desktop/Results/Lipidyzer")
R<-read.csv("Wake_Lipidyzer.csv",check.names=FALSE)
index<-apply(R[,1:10],1,is.na)
index<-apply(index,2,sum)
R<-R[index<1,]
index<- R$StudyID %in% c(1060202,1060127,1181043,1050103,1182023)
R<-R[index==FALSE,]
edata<-log2(t(R[,-c(1:25)]))
datExpr<-t(impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data)

powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- 16
adjacency <- adjacency(datExpr, power = softPower)
TOM <- TOMsimilarity(adjacency);
dissTOM <- 1-TOM
geneTree <- hclust(as.dist(dissTOM), method = "average");
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
minModuleSize <- 50
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
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
MEs_Lipidyzer <- MEList$eigengenes
colors<-MEList$validColors


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

powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- 9
adjacency <- adjacency(datExpr, power = softPower)
TOM <- TOMsimilarity(adjacency);
dissTOM <- 1-TOM
geneTree <- hclust(as.dist(dissTOM), method = "average");
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
minModuleSize <- 25
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
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
colors<-MEList$validColors


All_eigenfeatures<-cbind(as.data.frame(MEs_Lipidyzer),as.data.frame(MEs_Polar)) 
colnames(All_eigenfeatures)<-c("MEblack_L","MEblue_L","MEbrown_L","MEgreen_L","MEgrey_L","MEred_L","MEturquoise_L","MEyellow_L","MEgrey_P","MEmauve_P")

#######################Set Up Lipidyzer WGCNA 
setwd("~/Desktop/Results/Lipidyzer")
R<-read.csv("Wake_Lipidyzer.csv",check.names=FALSE)
index<-apply(R[,1:10],1,is.na)
index<-apply(index,2,sum)
R<-R[index<1,]
index<- R$StudyID %in% c(1060202,1060127,1181043,1050103,1182023)  
meta<-R[index==FALSE,c(2:11,14,18)]

final<-cbind(as.data.frame(meta),as.data.frame(All_eigenfeatures))


R<-read.csv("Gkk.csv",check.names=FALSE)[,-c(1:2)]
R$`Visit #`<-as.factor(R$`Visit #`)
R$Bsln_Cholesterol<-as.numeric(R$Bsln_Cholesterol)
R$Bsln_LDL<-as.numeric(R$Bsln_LDL)
R$Bsln_HDL<-as.numeric(R$Bsln_HDL)
R$Bsln_Trigly<-as.numeric(R$Bsln_Trigly)
R$Bsln_Fasting_Glu<-as.numeric(R$Bsln_Fasting_Glu)
T<-h2pc(R)
strength<-arc.strength(T,R)
strength.plot(T,strength)





combined<-cbind(as.data.frame(total),phate_clusters_energy,phate_clusters_lipid)
combined$`Visit #`<-as.factor(combined$`Visit #`)
combined$Bsln_Cholesterol<-as.numeric(combined$Bsln_Cholesterol)
combined$Bsln_LDL<-as.numeric(combined$Bsln_LDL)
combined$Bsln_HDL<-as.numeric(combined$Bsln_HDL)
combined$Bsln_Trigly<-as.numeric(combined$Bsln_Trigly)
combined$Bsln_Fasting_Glu<-as.numeric(combined$Bsln_Fasting_Glu)
combined$phate_clusters_lipid<-as.factor(combined$phate_clusters_lipid)
combined$phate_clusters_energy<-as.factor(combined$phate_clusters_energy)

qplot(data=combined,x=PHATE1,y=PHATE2,col=phate_clusters)
test<-combined%>%select(3:5,14,18,30:31)

combined1<-combined%>%filter(Study=="PACE 2" & `Visit #`=="1")%>%select(-c(1,2,12:13,15:17,19:29))
combined2<-combined%>%filter(Study=="PACE 2" & `Visit #`=="2")%>%select(-c(1,2,12:29))
combined3<-combined%>%filter(Study=="PACE-IGT" & `Visit #`=="1")%>%select(-c(1,2,12:29))
combined4<-combined%>%filter(Study=="PACE-IGT" & `Visit #`=="2")%>%select(-c(1,2,12:29))
combined5<-combined%>%filter(Study=="PACE-MCI" & `Visit #`=="1")%>%select(-c(1,2,12:29))
combined6<-combined%>%filter(Study=="PACE-MCI" & `Visit #`=="2")%>%select(-c(1,2,12:29))



net1<-h2pc(combined1)