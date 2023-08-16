library(dplyr)
library(impute)
library(limma)
library(Biobase)
library(sva)

####Read in LC-MS data
negative<-read.csv("mapstone_plasma_metabolomics_NEG.csv",check.names=FALSE,na.strings="0", header=FALSE)
rt_neg<-as.numeric(as.character(negative[2:dim(negative)[1],2]))
mz_neg<-as.numeric(as.character(negative[2:dim(negative)[1],1]))
neg_SampleID<-t(negative[1,3:dim(negative)[2]])
neg_abunds<-as.data.frame(t(negative[2:dim(negative)[1],3:dim(negative)[2]]))
colnames(neg_abunds)<-mz_neg
negativeIDs<-cbind(as.data.frame(neg_SampleID),as.data.frame(neg_abunds))
colnames(negativeIDs)[1]<-"Untargeted_Sample_ID"
negativeIDs$Untargeted_Sample_ID<-sub("01$","",negativeIDs$Untargeted_Sample_ID)

positive<-read.csv("mapstone_plasma_metabolomics_POS.csv",check.names=FALSE,na.strings="0", header=FALSE)
rt_pos<-as.numeric(as.character(positive[2:dim(positive)[1],2]))
mz_pos<-as.numeric(as.character(positive[2:dim(positive)[1],1]))
pos_SampleID<-t(positive[1,3:dim(positive)[2]])
pos_abunds<-as.data.frame(t(positive[2:dim(positive)[1],3:dim(positive)[2]]))
colnames(pos_abunds)<-mz_pos
positiveIDs<-cbind(as.data.frame(pos_SampleID),as.data.frame(pos_abunds))
colnames(positiveIDs)[1]<-"Untargeted_Sample_ID"
positiveIDs$Untargeted_Sample_ID<-sub("01$","",positiveIDs$Untargeted_Sample_ID)
colnames(positiveIDs)[duplicated(colnames(positiveIDs))==TRUE]<-lapply(colnames(positiveIDs)[duplicated(colnames(positiveIDs))==TRUE],paste,"1",sep=".")

####Bind to clinical data 
clinical<-read.csv("Complete_Wake_Pheno.csv",check.names=FALSE)[,-c(1)]%>%filter(Study=="PACE-MCI")
positive_final<-dplyr::inner_join(clinical,positiveIDs)
negative_final<-dplyr::inner_join(clinical,negativeIDs)
####Should evaluate TRUE 
all.equal(positive_final$StudyID,negative_final$StudyID)

####Impute abundances 
abunds<-apply(cbind(as.data.frame(positive_final[,-c(1:26)]),as.data.frame(negative_final[,-c(1:26)])),2,as.numeric)
colnames(abunds)<-seq(1,dim(abunds)[2])
T_abunds<-t(abunds)
edata<-impute.knn(T_abunds,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data

####Set up pheno object for ExpressionSet  
pheno_meta_prime<-as.data.frame(positive_final[,1:26])
index<-duplicated(pheno_meta_prime$StudyID)
pheno_meta<-pheno_meta_prime[index==FALSE,]
pheno_table<-pheno_meta[,c(1,2,13)]
colnames(pheno_table)<-c("StudyID","Treatment","Study")
Treat<-factor(paste(pheno_table$Study,pheno_table$Treatment,sep="."))
levels(Treat)<-list(#####"PACE2.Aerobic"=c("PACE 2.Aerobic"),
					#####"PACE2.Stretch"=c("PACE 2.Stretch"),
					#####"PACEIGT.Aerobic"=c("PACE-IGT.Aerobic"),
					#####"PACEIGT.Stretch"=c("PACE-IGT.Stretch"),
					"PACEMCI.Aerobic"=c("PACE-MCI.Aerobic"),
					"PACEMCI.Stretch"=c("PACE-MCI.Stretch"))
total<-cbind(Treat,pheno_table,pheno_meta[,3:10])
rownames(total)<-NULL
pheno_meta<-AnnotatedDataFrame(total)

####Collapse pre-post variance into log-normalized change score
all.equal(pheno_meta_prime$StudyID[pheno_meta_prime$`Visit #`=="1"],
	      pheno_meta_prime$StudyID[pheno_meta_prime$`Visit #`=="2"])
time0<-edata[,pheno_meta_prime$`Visit #`=="1"]
time1<-edata[,pheno_meta_prime$`Visit #`=="2"]
squared<-(time1-time0)^2
edata<-sqrt(squared/(time1+time0))
rownames(edata)<-rownames(T_abunds)
colnames(edata)<-rownames(total)
edata<-log2(edata)

####Set up feature data for ExpressionSet 
mode<-as.factor(c(rep("ESI+",dim(positive_final[,-c(1:26)])[2]),rep("ESI-",dim(negative_final[,-c(1:26)])[2])))
feature_meta<-cbind(rbind(cbind(mz_pos,rt_pos),cbind(mz_neg,rt_neg)),as.data.frame(as.factor(mode)))
rownames(feature_meta)<-rownames(T_abunds)
colnames(feature_meta)<-c("MZ","RT","Mode")

####Features for Output 
feature_meta_use<-cbind(rownames(T_abunds),feature_meta)
colnames(feature_meta_use)[1]<-"Feature"

####Write feature object 
feature_meta<-AnnotatedDataFrame(feature_meta)

####Construct final ExpressionData object
combined<-ExpressionSet(assayData=edata,phenoData=pheno_meta,featureData=feature_meta)

####Fit SVs 
mod<-model.matrix(~as.factor(Treat),data=combined)
mod0<-model.matrix(~1,data=combined)
n.sv_LEEK<-num.sv(edata,mod,seed=122,method="leek")
n.sv_BE<-num.sv(edata,mod,seed=122)
svobj_LEEK<-sva(edata,mod,mod0,n.sv=n.sv_LEEK)$sv
svobj_BE<-sva(edata,mod,mod0,n.sv=n.sv_BE)$sv

####Write MCI SVs (BE+LEEK) to CSV
A<-cbind(total[,1:2],svobj_BE)
colnames(A)[3:6]<-c("SV1.BE","SV2.BE","SV3.BE","SV4.BE")
write.csv(A,file="PACE_MCI_SVs_METHOD_BE.csv")

B<-cbind(total[,1:2],svobj_LEEK)
colnames(B)[3:5]<-c("SV1.LEEK","SV2.LEEK","SV3.LEEK")
write.csv(B,file="PACE_MCI_SVs_METHOD_LEEK.csv")