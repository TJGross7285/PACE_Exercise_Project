library(dplyr)
library(impute)
library(limma)
library(Biobase)
library(sva)

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
####squared<-(time1-time0)^2
#####edata<-sqrt(squared/(time1+time0))
edata<-(time1-time0)/(time1+time0)
rownames(edata)<-rownames(T_abunds)
colnames(edata)<-rownames(total)
#####edata<-log2(edata)

mode<-as.factor(c(rep("ESI+",dim(positive_final[,-c(1:26)])[2]),rep("ESI-",dim(negative_final[,-c(1:26)])[2])))
feature_meta<-cbind(rbind(cbind(mz_pos,rt_pos),cbind(mz_neg,rt_neg)),as.data.frame(as.factor(mode)))
rownames(feature_meta)<-rownames(T_abunds)
colnames(feature_meta)<-c("MZ","RT","Mode")

####Features for Output 
feature_meta_use<-cbind(rownames(T_abunds),feature_meta)
colnames(feature_meta_use)[1]<-"Feature"

####Fit SVs BE 
mod<-model.matrix(~as.factor(Treat),data=total)
mod0<-model.matrix(~1,data=total)
n.sv_BE<-num.sv(edata,mod,seed=122)
svobj_BE<-sva(edata,mod,mod0,n.sv=n.sv_BE)$sv


modSv<-cbind(mod,svobj_BE)
mod0Sv<-cbind(mod0,svobj_BE)
pValuesSv<-f.pvalue(edata,modSv,mod0Sv)

fc<-apply(edata[,total$Treat=="PACEMCI.Aerobic"],1,median)/apply(edata[,total$Treat=="PACEMCI.Stretch"],1,median)

table<-cbind(feature_meta_use,pValuesSv,fc)
colnames(table)[c(5,6)]<-c("P.Value","FC_Aerobic.Stretch")

write.table(table%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","FC_Aerobic.Stretch"),
			sep="\t",file="PACE_5_20_2021.POS.txt",row.names=FALSE)
write.table(table%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","FC_Aerobic.Stretch"),
			sep="\t",file="PACE_5_20_2021.NEG.txt",row.names=FALSE)



cd /Volumes/NO\ NAME/Active_Projects_4_9_2021/PACE_4_2021 

mummichog -f PACE_5_20_2021.NEG.txt -o PACE_5_20_2021.NEG -m negative
mummichog -f PACE_5_20_2021.POS.txt -o PACE_5_20_2021.POS -m positive
 












############ORIGINAL 5/20/2021
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

####Fit SVs BE 
mod<-model.matrix(~as.factor(Treat),data=combined)
mod0<-model.matrix(~1,data=combined)
n.sv_BE<-num.sv(edata,mod,seed=122)
svobj_BE<-sva(edata,mod,mod0,n.sv=n.sv_BE)$sv

design_table<-cbind(model.matrix(~0+as.factor(Treat),data=combined),as.data.frame(svobj_BE)) #########svobj_BE$sv

colnames(design_table)[1:2]<-c("PACEMCI.Aerobic","PACEMCI.Stretch")   ############"PACE2.Aerobic","PACE2.Stretch","PACEIGT.Aerobic","PACEIGT.Stretch",
						  
arrayw<-arrayWeights(edata, design=design_table)

####Fit linear model
fit1<-lmFit(edata,design_table, weights=arrayw)

cm3<-makeContrasts(
	`PACEMCI.Aerobic-PACEMCI.Stretch`= PACEMCI.Aerobic-PACEMCI.Stretch, 
	levels=design_table)


fit3_F <- contrasts.fit(fit1, cm3)
fit3_F <- eBayes(fit3_F,trend=TRUE)
T<-topTableF(fit3_F,adjust="BH",number=100000)
V<-cbind(rownames(T),T)
colnames(V)[1]<-"Feature"
joinV<-dplyr::inner_join(feature_meta_use,V)
write.table(joinV%>%filter(Mode=="ESI+")%>%dplyr::select("MZ","RT","P.Value","PACEMCI.Aerobic.PACEMCI.Stretch"),
			sep="\t",file="PACE_MCI.POS_BE.txt",row.names=FALSE)
write.table(joinV%>%filter(Mode=="ESI-")%>%dplyr::select("MZ","RT","P.Value","PACEMCI.Aerobic.PACEMCI.Stretch"),
			sep="\t",file="PACE_MCI.NEG_BE.txt",row.names=FALSE)












#########################################################
###################Fit SVs Leek
mod<-model.matrix(~as.factor(Treat),data=combined)
mod0<-model.matrix(~1,data=combined)
n.sv_LEEK<-num.sv(edata,mod,seed=122, method="leek")
svobj_LEEK<-sva(edata,mod,mod0,n.sv=n.sv_LEEK)$sv

design_table<-cbind(model.matrix(~0+as.factor(Treat),data=combined),as.data.frame(svobj_LEEK))

colnames(design_table)[1:2]<-c("PACEMCI.Aerobic","PACEMCI.Stretch")   ############"PACE2.Aerobic","PACE2.Stretch","PACEIGT.Aerobic","PACEIGT.Stretch",
						  
arrayw<-arrayWeights(edata, design=design_table)

####Fit linear model
fit1<-lmFit(edata,design_table, weights=arrayw)

cm3<-makeContrasts(
	`PACEMCI.Aerobic-PACEMCI.Stretch`= PACEMCI.Aerobic-PACEMCI.Stretch, 
	levels=design_table)


fit3_F <- contrasts.fit(fit1, cm3)
fit3_F <- eBayes(fit3_F,trend=TRUE)
T<-topTableF(fit3_F,adjust="BH",number=100000)
V<-cbind(rownames(T),T)
colnames(V)[1]<-"Feature"
joinZ<-dplyr::inner_join(feature_meta_use,V)
write.table(joinZ%>%filter(Mode=="ESI+")%>%dplyr::select("MZ","RT","P.Value","PACEMCI.Aerobic.PACEMCI.Stretch"),
			sep="\t",file="PACE_MCI.POS_LEEK.txt",row.names=FALSE)
write.table(joinZ%>%filter(Mode=="ESI-")%>%dplyr::select("MZ","RT","P.Value","PACEMCI.Aerobic.PACEMCI.Stretch"),
			sep="\t",file="PACE_MCI.NEG_LEEK.txt",row.names=FALSE)

###########Bash Mummichog
cd /Volumes/NO\ NAME/Active_Projects_4_9_2021/PACE_4_2021 

mummichog -f PACE_MCI.NEG_BE.txt -o PACE_MCI.NEG_BE -m negative
mummichog -f PACE_MCI.POS_BE.txt -o PACE_MCI.POS_BE -m positive

mummichog -f PACE_MCI.NEG_LEEK.txt -o PACE_MCI.NEG_LEEK -m negative
mummichog -f PACE_MCI.POS_LEEK.txt -o PACE_MCI.POS_LEEK -m positive


#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################

library(dplyr)
library(impute)
library(limma)
library(Biobase)
library(sva)

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
clinical<-read.csv("Complete_Wake_Pheno.csv",check.names=FALSE)[,-c(1)]%>%filter(Study=="PACE-IGT")
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
levels(Treat)<-list("PACE2.Aerobic"=c("PACE 2.Aerobic"),
					"PACE2.Stretch"=c("PACE 2.Stretch"),
					"PACEIGT.Aerobic"=c("PACE-IGT.Aerobic"),
					"PACEIGT.Stretch"=c("PACE-IGT.Stretch"),
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

####Fit SVs BE 
mod<-model.matrix(~as.factor(factor(Treat)),data=combined)
mod0<-model.matrix(~1,data=combined)
n.sv_BE<-num.sv(edata,mod,seed=122)
svobj_BE<-sva(edata,mod,mod0,n.sv=n.sv_BE)$sv

design_table<-cbind(model.matrix(~0+as.factor(factor(Treat)),data=combined),as.data.frame(svobj_BE)) #########svobj_BE$sv

colnames(design_table)[1:2]<-c("PACEIGT.Aerobic","PACEIGT.Stretch")
	                           
						  
arrayw<-arrayWeights(edata, design=design_table)

####Fit linear model
fit1<-lmFit(edata,design_table, weights=arrayw)    

cm2<-makeContrasts(`PACEIGT.Aerobic-PACEIGT.Stretch` = PACEIGT.Aerobic-PACEIGT.Stretch,
	               levels=design_table)


fit2_F <- contrasts.fit(fit1, cm2)
fit2_F <- eBayes(fit2_F,trend=TRUE)
T<-topTableF(fit2_F,adjust="BH",number=100000)
U<-cbind(rownames(T),T)
colnames(U)[1]<-"Feature"
joinU<-dplyr::inner_join(feature_meta_use,U)
write.table(joinU%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","PACEIGT.Aerobic.PACEIGT.Stretch"),
			sep="\t",file="PACE_IGT.POS_BE.txt",row.names=FALSE)
write.table(joinU%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","PACEIGT.Aerobic.PACEIGT.Stretch"),
			sep="\t",file="PACE_IGT.NEG_BE.txt",row.names=FALSE)








#########################################################
###################Fit SVs Leek
mod<-model.matrix(~as.factor(factor(Treat)),data=combined)
mod0<-model.matrix(~1,data=combined)
n.sv_LEEK<-num.sv(edata,mod,seed=122, method="leek")
svobj_LEEK<-sva(edata,mod,mod0,n.sv=n.sv_LEEK)$sv

design_table<-cbind(model.matrix(~0+as.factor(factor(Treat)),data=combined),as.data.frame(svobj_LEEK))

colnames(design_table)[1:2]<-c("PACEIGT.Aerobic","PACEIGT.Stretch")
						  
arrayw<-arrayWeights(edata, design=design_table)

####Fit linear model
fit1<-lmFit(edata,design_table, weights=arrayw)

cm2<-makeContrasts(
	`PACEIGT.Aerobic-PACEIGT.Stretch` = PACEIGT.Aerobic-PACEIGT.Stretch, 
	levels=design_table)


fit2_F <- contrasts.fit(fit1, cm2)
fit2_F <- eBayes(fit2_F,trend=TRUE)
T<-topTableF(fit2_F,adjust="BH",number=100000)
U<-cbind(rownames(T),T)
colnames(U)[1]<-"Feature"
joinV<-dplyr::inner_join(feature_meta_use,U)
write.table(joinV%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","PACEIGT.Aerobic.PACEIGT.Stretch"),
			sep="\t",file="PACE_IGT.POS_LEEK.txt",row.names=FALSE)
write.table(joinV%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","PACEIGT.Aerobic.PACEIGT.Stretch"),
			sep="\t",file="PACE_IGT.NEG_LEEK.txt",row.names=FALSE)

###########Bash Mummichog
cd /Volumes/NO\ NAME/Active_Projects_4_9_2021/PACE_4_2021 

mummichog -f PACE_IGT.NEG_BE.txt -o PACE_IGT.NEG_BE -m negative
mummichog -f PACE_IGT.POS_BE.txt -o PACE_IGT.POS_BE -m positive

mummichog -f PACE_IGT.NEG_LEEK.txt -o PACE_IGT.NEG_LEEK -m negative
mummichog -f PACE_IGT.POS_LEEK.txt -o PACE_IGT.POS_LEEK -m positive




#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################

library(dplyr)
library(impute)
library(limma)
library(Biobase)
library(sva)

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
clinical<-read.csv("Complete_Wake_Pheno.csv",check.names=FALSE)[,-c(1)]%>%filter(Study=="PACE 2")
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
levels(Treat)<-list("PACE2.Aerobic"=c("PACE 2.Aerobic"),
					"PACE2.Stretch"=c("PACE 2.Stretch"),
					"PACEIGT.Aerobic"=c("PACE-IGT.Aerobic"),
					"PACEIGT.Stretch"=c("PACE-IGT.Stretch"),
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

####Fit SVs BE 
mod<-model.matrix(~as.factor(factor(Treat)),data=combined)
mod0<-model.matrix(~1,data=combined)
n.sv_BE<-num.sv(edata,mod,seed=122)
svobj_BE<-sva(edata,mod,mod0,n.sv=n.sv_BE)$sv

design_table<-cbind(model.matrix(~0+as.factor(factor(Treat)),data=combined),as.data.frame(svobj_BE)) #########svobj_BE$sv

colnames(design_table)[1:2]<-c("PACE2.Aerobic","PACE2.Stretch")
	                           
						  
arrayw<-arrayWeights(edata, design=design_table)

####Fit linear model
fit1<-lmFit(edata,design_table, weights=arrayw)    

cm2<-makeContrasts(`PACE2.Aerobic-PACE2.Stretch` = PACE2.Aerobic-PACE2.Stretch,
	               levels=design_table)


fit2_F <- contrasts.fit(fit1, cm2)
fit2_F <- eBayes(fit2_F,trend=TRUE)
T<-topTableF(fit2_F,adjust="BH",number=100000)
U<-cbind(rownames(T),T)
colnames(U)[1]<-"Feature"
joinU<-dplyr::inner_join(feature_meta_use,U)
write.table(joinU%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","PACE2.Aerobic.PACE2.Stretch"),
			sep="\t",file="PACE_2.POS_BE.txt",row.names=FALSE)
write.table(joinU%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","PACE2.Aerobic.PACE2.Stretch"),
			sep="\t",file="PACE_2.NEG_BE.txt",row.names=FALSE)








#########################################################
###################Fit SVs Leek
mod<-model.matrix(~as.factor(factor(Treat)),data=combined)
mod0<-model.matrix(~1,data=combined)
n.sv_LEEK<-num.sv(edata,mod,seed=122, method="leek")
svobj_LEEK<-sva(edata,mod,mod0,n.sv=n.sv_LEEK)$sv

design_table<-cbind(model.matrix(~0+as.factor(factor(Treat)),data=combined),as.data.frame(svobj_LEEK))

colnames(design_table)[1:2]<-c("PACE2.Aerobic","PACE2.Stretch")
						  
arrayw<-arrayWeights(edata, design=design_table)

####Fit linear model
fit1<-lmFit(edata,design_table, weights=arrayw)

cm2<-makeContrasts(
	`PACE2.Aerobic-PACE2.Stretch` = PACE2.Aerobic-PACE2.Stretch, 
	levels=design_table)


fit2_F <- contrasts.fit(fit1, cm2)
fit2_F <- eBayes(fit2_F,trend=TRUE)
T<-topTableF(fit2_F,adjust="BH",number=100000)
U<-cbind(rownames(T),T)
colnames(U)[1]<-"Feature"
joinV<-dplyr::inner_join(feature_meta_use,U)
write.table(joinV%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","PACE2.Aerobic.PACE2.Stretch"),
			sep="\t",file="PACE_2.POS_LEEK.txt",row.names=FALSE)
write.table(joinV%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","PACE2.Aerobic.PACE2.Stretch"),
			sep="\t",file="PACE_2.NEG_LEEK.txt",row.names=FALSE)

###########Bash Mummichog
cd /Volumes/NO\ NAME/Active_Projects_4_9_2021/PACE_4_2021 

mummichog -f PACE_2.NEG_BE.txt -o PACE_2.NEG_BE -m negative
mummichog -f PACE_2.POS_BE.txt -o PACE_2.POS_BE -m positive

mummichog -f PACE_2.NEG_LEEK.txt -o PACE_2.NEG_LEEK -m negative
mummichog -f PACE_2.POS_LEEK.txt -o PACE_2.POS_LEEK -m positive


mummichog -f PACE_IGT.NEG_BE.txt -o PACE_IGT.NEG_BE -m negative
mummichog -f PACE_IGT.POS_BE.txt -o PACE_IGT.POS_BE -m positive

mummichog -f PACE_IGT.NEG_LEEK.txt -o PACE_IGT.NEG_LEEK -m negative
mummichog -f PACE_IGT.POS_LEEK.txt -o PACE_IGT.POS_LEEK -m positive


mummichog -f PACE_MCI.NEG_BE.txt -o PACE_MCI.NEG_BE -m negative
mummichog -f PACE_MCI.POS_BE.txt -o PACE_MCI.POS_BE -m positive

mummichog -f PACE_MCI.NEG_LEEK.txt -o PACE_MCI.NEG_LEEK -m negative
mummichog -f PACE_MCI.POS_LEEK.txt -o PACE_MCI.POS_LEEK -m positive




