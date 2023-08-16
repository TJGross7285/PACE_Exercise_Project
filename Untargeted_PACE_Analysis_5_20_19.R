

library(dplyr)
library(impute)
library(limma)

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
U<-read.csv("mapstone_plasma_metabolomics_NEG.csv",check.names=FALSE)
I<-read.csv("mapstone_plasma_metabolomics_POS.csv",check.names=FALSE)

####Format negative mode data
abunds_neg<-cbind(as.data.frame(colnames(U)[-c(1:2)]),as.data.frame(t(U[,-c(1:2)])))
colnames(abunds_neg)[1]<-"Untargeted_Sample_ID"
abunds_neg$Untargeted_Sample_ID<-sub("01$","",abunds_neg$Untargeted_Sample_ID)
join_neg<-dplyr::inner_join(join,abunds_neg,by="Untargeted_Sample_ID")
colnames(join_neg)[18:length(join_neg)]<-paste0(colnames(join_neg)[18:length(join_neg)],rep(".NEG",length(colnames(join_neg)[18:length(join_neg)])))

####Format postive mode data
abunds_pos<-cbind(as.data.frame(colnames(I)[-c(1:2)]),as.data.frame(t(I[,-c(1:2)])))
colnames(abunds_pos)[1]<-"Untargeted_Sample_ID"
abunds_pos$Untargeted_Sample_ID<-sub("01$","",abunds_pos$Untargeted_Sample_ID)
join_pos<-dplyr::inner_join(join,abunds_pos,by="Untargeted_Sample_ID")
colnames(join_pos)[18:length(join_pos)]<-paste0(colnames(join_pos)[18:length(join_pos)],rep(".POS",length(colnames(join_pos)[18:length(join_pos)])))

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

####Set up accessory objects for DE incorporating pre-post timepoints
pheno_table<-total_frame[,c(1,2,13,17)]
colnames(pheno_table)<-c("StudyID","Treatment","Study","Main")
Treat<-factor(paste(pheno_table$Study,pheno_table$Main,sep="."))
Treat<-factor(paste(Treat,pheno_table$Treatment,sep="."))
levels(Treat)<-list("PACE2pre.Aerobic"=c("PACE 2.1.Aerobic"),
					"PACE2post.Aerobic"=c("PACE 2.2.Aerobic"),
					"PACE2pre.Stretching"=c("PACE 2.1.Stretch"),
					"PACE2post.Stretching"=c("PACE 2.2.Stretch"),
					"PACEIGTpre.Aerobic"=c("PACE-IGT.1.Aerobic"),
					"PACEIGTpost.Aerobic"=c("PACE-IGT.2.Aerobic"),
					"PACEIGTpre.Stretching"=c("PACE-IGT.1.Stretch"),
					"PACEIGTpost.Stretching"=c("PACE-IGT.2.Stretch"),
					"PACEMCIpre.Aerobic"=c("PACE-MCI.1.Aerobic"),
					"PACEMCIpost.Aerobic"=c("PACE-MCI.2.Aerobic"),
					"PACEMCIpre.Stretching"=c("PACE-MCI.1.Stretch"),
					"PACEMCIpost.Stretching"=c("PACE-MCI.2.Stretch"))
total<-cbind(Treat,pheno_table,total_frame[,3:10])
####+Sex+APOE+Bsln_BMI+Bsln_Cholesterol+Bsln_LDL+Bsln_HDL+Bsln_Trigly+Bsln_Fasting_Glu
design1<-model.matrix(~0+Treat,data=total)
colnames(design1)[1:12]<-c("PACE2pre.Aerobic",
						"PACE2post.Aerobic",
						"PACE2pre.Stretching",
						"PACE2post.Stretching",
						"PACEIGTpre.Aerobic",
						"PACEIGTpost.Aerobic",
						"PACEIGTpre.Stretching",
						"PACEIGTpost.Stretching",
						"PACEMCIpre.Aerobic",
						"PACEMCIpost.Aerobic",
						"PACEMCIpre.Stretching",
						"PACEMCIpost.Stretching")
corfit1<-duplicateCorrelation(edata,design1,block=pheno_table$StudyID)

arrayw<-arrayWeights(edata, design=design1)

####Fit linear model for within subjects contrast
fit1<-lmFit(edata,design1,block=pheno_table$StudyID,correlation=corfit1$consensus,weights=arrayw)
cm<- makeContrasts(
	`PACE2post.Aerobic-PACE2pre.Aerobic` = PACE2post.Aerobic-PACE2pre.Aerobic,
	`PACEIGTpost.Aerobic-PACEIGTpre.Aerobic` = PACEIGTpost.Aerobic-PACEIGTpre.Aerobic,
	`PACEMCIpost.Aerobic-PACEMCIpre.Aerobic` = PACEMCIpost.Aerobic-PACEMCIpre.Aerobic,
	`PACE2post.Stretching-PACE2pre.Stretching` = PACE2post.Stretching-PACE2pre.Stretching,
	`PACEIGTpost.Stretching-PACEIGTpre.Stretching` = PACEIGTpost.Stretching-PACEIGTpre.Stretching,
	`PACEMCIpost.Stretching-PACEMCIpre.Stretching` = PACEMCIpost.Stretching-PACEMCIpre.Stretching,
	levels=design1)

cm1<-makeContrasts(`PACE2post.Aerobic-PACE2pre.Aerobic` = PACE2post.Aerobic-PACE2pre.Aerobic,levels=design1)
cm2<-makeContrasts(`PACEIGTpost.Aerobic-PACEIGTpre.Aerobic` = PACEIGTpost.Aerobic-PACEIGTpre.Aerobic,levels=design1)
cm3<-makeContrasts(`PACEMCIpost.Aerobic-PACEMCIpre.Aerobic` = PACEMCIpost.Aerobic-PACEMCIpre.Aerobic,levels=design1)

cm4<-makeContrasts(`PACE2post.Stretching-PACE2pre.Stretching` = PACE2post.Stretching-PACE2pre.Stretching,levels=design1)
cm5<-makeContrasts(`PACEIGTpost.Stretching-PACEIGTpre.Stretching` = PACEIGTpost.Stretching-PACEIGTpre.Stretching,levels=design1)
cm6<-makeContrasts(`PACEMCIpost.Stretching-PACEMCIpre.Stretching` = PACEMCIpost.Stretching-PACEMCIpre.Stretching,levels=design1)

fit1_F <- contrasts.fit(fit1, cm1)
fit1_F <- eBayes(fit1_F,trend=TRUE)

fit2_F<- contrasts.fit(fit1, cm2)
fit2_F<-eBayes(fit2_F,trend=TRUE)

fit3_F<-contrasts.fit(fit1, cm3)
fit3_F<-eBayes(fit3_F,trend=TRUE)


T<-topTableF(fit3_F,number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"FeatureID"
Tplus<-dplyr::inner_join(featureData,T)%>%filter(Mode=="ESI+")%>%dplyr::select("MZ","RT","P.Value","F")
Tminus<-dplyr::inner_join(featureData,T)%>%filter(Mode=="ESI-")%>%dplyr::select("MZ","RT","P.Value","F")

write.table(Tplus, sep="\t",file="PACE2_POS.txt",row.names=FALSE)
write.table(Tminus, sep="\t",file="PACE2_NEG.txt",row.names=FALSE)


T<-topTableF(fit2_F,number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"FeatureID"
Tplus<-dplyr::inner_join(featureData,T)%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","F")
Tminus<-dplyr::inner_join(featureData,T)%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","F")

write.table(Tplus, sep="\t",file="PACEMCI_POS.txt",row.names=FALSE)
write.table(Tminus, sep="\t",file="PACEMCI_NEG.txt",row.names=FALSE)



fit1_F <- contrasts.fit(fit1, cm1)
fit1_F <- eBayes(fit1_F,trend=TRUE)

fit2_F<- contrasts.fit(fit1, cm2)
fit2_F<-eBayes(fit2_F,trend=TRUE)

fit3_F<-contrasts.fit(fit1, cm3)
fit3_F<-eBayes(fit3_F,trend=TRUE)

fit4_F <- contrasts.fit(fit1, cm4)
fit4_F <- eBayes(fit4_F,trend=TRUE)

fit5_F<- contrasts.fit(fit1, cm5)
fit5_F<-eBayes(fit5_F,trend=TRUE)

fit6_F<-contrasts.fit(fit1, cm6)
fit6_F<-eBayes(fit6_F,trend=TRUE)


A<-topTableF(fit1_F,number=100000)
T<-cbind(rownames(A),A)
colnames(T)[1]<-"FeatureID"
A<-dplyr::inner_join(featureData,T)%>%arrange(P.Value)

B<-topTableF(fit2_F,number=100000)
T<-cbind(rownames(B),B)
colnames(T)[1]<-"FeatureID"
B<-dplyr::inner_join(featureData,T)%>%arrange(P.Value)

C<-topTableF(fit3_F,number=100000)
T<-cbind(rownames(C),C)
colnames(T)[1]<-"FeatureID"
C<-dplyr::inner_join(featureData,T)%>%arrange(P.Value)

D<-topTableF(fit4_F,number=100000)
T<-cbind(rownames(D),D)
colnames(T)[1]<-"FeatureID"
D<-dplyr::inner_join(featureData,T)%>%arrange(P.Value)

E<-topTableF(fit5_F,number=100000)
T<-cbind(rownames(E),E)
colnames(T)[1]<-"FeatureID"
E<-dplyr::inner_join(featureData,T)%>%arrange(P.Value)

F<-topTableF(fit6_F,number=100000)
T<-cbind(rownames(F),F)
colnames(T)[1]<-"FeatureID"
F<-dplyr::inner_join(featureData,T)%>%arrange(P.Value)


pdf(file="PACE_Volcano_12_23.pdf")
plot(A$`PACE2post.Aerobic.PACE2pre.Aerobic`,-log10(A$P.Value))
abline(a=1.30103,b=0,col="firebrick3")
plot(B$`PACEIGTpost.Aerobic.PACEIGTpre.Aerobic`,-log10(B$P.Value))
abline(a=1.30103,b=0,col="firebrick3")
plot(C$`PACEMCIpost.Aerobic.PACEMCIpre.Aerobic`,-log10(C$P.Value))
abline(a=1.30103,b=0,col="firebrick3")
plot(D$`PACE2post.Stretching.PACE2pre.Stretching`,-log10(D$P.Value))
abline(a=1.30103,b=0,col="firebrick3")
plot(E$`PACEIGTpost.Stretching.PACEIGTpre.Stretching`,-log10(E$P.Value))
abline(a=1.30103,b=0,col="firebrick3")
plot(F$`PACEMCIpost.Stretching.PACEMCIpre.Stretching`,-log10(F$P.Value))
abline(a=1.30103,b=0,col="firebrick3")
dev.off()
















U<-topTableF(fit2_F,number=100000)
V<-topTableF(fit3_F,number=100000)









Y<-topTableF(fit1_F,coef=2,number=100000)
Y_anno<-rep("PACEIGTpost.Aerobic-PACEIGTpre.Aerobic",dim(Y)[1])
correctedFC<-predFCm(fit1_F,coef=2)
Y<-cbind(Y,correctedFC)

U<-topTable(fit1_F,coef=3,number=100000)
U_anno<-rep("PACEMCIpost.Aerobic-PACEMCIpre.Aerobic",dim(U)[1])
correctedFC<-predFCm(fit1_F,coef=3)
U<-cbind(U,correctedFC)

I<-topTableF(fit1_F,coef=4,number=100000)
I_anno<-rep("PACE2post.Stretching-PACE2pre.Stretching",dim(I)[1])
correctedFC<-predFCm(fit1_F,coef=4)
I<-cbind(I,correctedFC)

O<-topTableF(fit1_F,coef=5,number=100000)
O_anno<-rep("PACEIGTpost.Stretching-PACEIGTpre.Stretching",dim(O)[1])
correctedFC<-predFCm(fit1_F,coef=5)
O<-cbind(O,correctedFC)

P<-topTableF(fit1_F,coef=6,number=100000)
P_anno<-rep("PACEMCIpost.Stretching-PACEMCIpre.Stretching",dim(P)[1])
correctedFC<-predFCm(fit1_F,coef=6)
P<-cbind(P,correctedFC)

study_vector<-c(T_anno,Y_anno,U_anno,I_anno,O_anno,P_anno)
feature_vector<-c(rownames(T),rownames(Y),rownames(U),rownames(I),rownames(O),rownames(P))
total_tables<-rbind(T,Y,U,I,O,P)
final_frame<-cbind(feature_vector,study_vector,total_tables)
colnames(final_frame)[1]<-"Feature"
colnames(final_frame)[2]<-"Contrast"

write.csv(final_frame, file="PACE_Untargeted_6_12_PRE_POST_Contrasts.csv")



list1<-final_frame%>%filter(Contrast=="PACE2post.Aerobic-PACE2pre.Aerobic" & P.Value<.05 & abs(logFC)>.1)%>%select(Feature,logFC)
list2<-final_frame%>%filter(Contrast=="PACEIGTpost.Aerobic-PACEIGTpre.Aerobic" & P.Value<.05 & abs(logFC)>.1)%>%select(Feature,logFC)
list3<-final_frame%>%filter(Contrast=="PACEMCIpost.Aerobic-PACEMCIpre.Aerobic" & P.Value<.05 & abs(logFC)>.1)%>%select(Feature,logFC)
list4<-final_frame%>%filter(Contrast=="PACE2post.Stretching-PACE2pre.Stretching" & P.Value<.05 & abs(logFC)>.1)%>%select(Feature,logFC)
list5<-final_frame%>%filter(Contrast=="PACEIGTpost.Stretching-PACEIGTpre.Stretching" & P.Value<.05 & abs(logFC)>.1)%>%select(Feature,logFC)
list6<-final_frame%>%filter(Contrast=="PACEMCIpost.Stretching-PACEMCIpre.Stretching" & P.Value<.05 & abs(logFC)>.1)%>%select(Feature,logFC)


#********NONE*********
intersect(list1$Feature,list4$Feature)

#********NONE*********
intersect(list2$Feature,list5$Feature)

#********HITS*********
intersect(list3$Feature,list6$Feature)

pdf(file="DE_Untargeted_PACE_Venn_PRE_POST.pdf")

GOVenn(list1,list2,list3,label=c("PACE2post.Aerobic-PACE2pre.Aerobic",
								"PACEIGTpost.Aerobic-PACEIGTpre.Aerobic",
								"PACEMCIpost.Aerobic-PACEMCIpre.Aerobic"),
								lfc.col=c("seagreen","gold","firebrick2"))

GOVenn(list4,list5,list6,label=c("PACE2post.Stretching-PACE2pre.Stretching",
								"PACEIGTpost.Stretching-PACEIGTpre.Stretching",
								"PACEMCIpost.Stretching-PACEMCIpre.Stretching"),
								lfc.col=c("seagreen","gold","firebrick2"))

GOVenn(list1,list4,label=c("PACE2post.Aerobic-PACE2pre.Aerobic",
								"PACE2post.Stretching-PACE2pre.Stretching"),
								lfc.col=c("seagreen","gold","firebrick2"))

GOVenn(list2,list5,label=c("PACEIGTpost.Aerobic-PACEIGTpre.Aerobic",
								"PACEIGTpost.Stretching-PACEIGTpre.Stretching"),
								lfc.col=c("seagreen","gold","firebrick2"))

GOVenn(list3,list6,label=c("PACEMCIpost.Aerobic-PACEMCIpre.Aerobic",
							"PACEMCIpost.Stretching-PACEMCIpre.Stretching"),
								lfc.col=c("seagreen","gold","firebrick2"))
dev.off()



























####Set up accessory objects for DE
pheno_table<-join[,c(1,2,13,17)]
colnames(pheno_table)<-c("StudyID","Treatment","Study","Main")
Treat<-factor(paste(pheno_table$Study,pheno_table$Main,sep="."))
Treat<-factor(paste(Treat,pheno_table$Treatment,sep="."))
levels(Treat)<-list("PACE2pre.Aerobic"=c("PACE 2.1.Aerobic"),
					"PACE2post.Aerobic"=c("PACE 2.2.Aerobic"),
					"PACE2pre.Stretch"=c("PACE 2.1.Stretch"),
					"PACE2post.Stretch"=c("PACE 2.2.Stretch"),
					"PACEIGTpre.Aerobic"=c("PACE-IGT.1.Aerobic"),
					"PACEIGTpost.Aerobic"=c("PACE-IGT.2.Aerobic"),
					"PACEIGTpre.Stretch"=c("PACE-IGT.1.Stretch"),
					"PACEIGTpost.Stretch"=c("PACE-IGT.2.Stretch"),
					"PACEMCIpre.Aerobic"=c("PACE-MCI.1.Aerobic"),
					"PACEMCIpost.Aerobic"=c("PACE-MCI.2.Aerobic"),
					"PACEMCIpre.Stretch"=c("PACE-MCI.1.Stretch"),
					"PACEMCIpost.Stretch"=c("PACE-MCI.2.Stretch"))
design<-model.matrix(~0+Treat)
colnames(design)<-levels(Treat)
corfit<-duplicateCorrelation(edata,design,block=pheno_table$StudyID)

####Fit linear model for within subjects contrast
fit<-lmFit(edata,design,block=pheno_table$StudyID,correlation=corfit$consensus)
cm <- makeContrasts(
	`PACE2post.Aerobic-PACE2pre.Aerobic` = PACE2post.Aerobic-PACE2pre.Aerobic,
	`PACEIGTpost.Aerobic-PACEIGTpre.Aerobic` = PACEIGTpost.Aerobic-PACEIGTpre.Aerobic,
	`PACEMCIpost.Aerobic-PACEMCIpre.Aerobic` = PACEMCIpost.Aerobic-PACEMCIpre.Aerobic,
	`PACE2post.Stretch-PACE2pre.Stretch` = PACE2post.Stretch-PACE2pre.Stretch,
	`PACEIGTpost.Stretch-PACEIGTpre.Stretch` = PACEIGTpost.Stretch-PACEIGTpre.Stretch,
	`PACEMCIpost.Stretch-PACEMCIpre.Stretch` = PACEMCIpost.Stretch-PACEMCIpre.Stretch,
	levels=design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
T<-topTable(fit2,number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"FeatureID"

image(sub)



####Write DE table to file
Final_DE<-dplyr::inner_join(as.data.frame(featureData),T,by="FeatureID")%>%filter(Mode=="-")%>%select("MZ","RT","P.Value","F")
write.csv(Final_DE,file="Wake_Final_DE_Untargeted_WITHIN_NEG.csv")
Final_DE<-dplyr::inner_join(as.data.frame(featureData),T,by="FeatureID")%>%filter(Mode=="+")%>%select("MZ","RT","P.Value","F")
write.csv(Final_DE,file="Wake_Final_DE_Untargeted_WITHIN_POS.csv")







####Calculate correlations of metabolites to other pheno info 
library(gplots)
pheno_corr<-apply(total_frame[,c(5:10)],2,as.numeric)
abunds<-total_frame[,-c(1:26)]
corr<-cor(abunds,pheno_corr,method="spearman",use="complete.obs")
pdf(file="Untargeted_Pheno_Corr.pdf")
heatmap.2(corr,trace="none")
dev.off()












data_lip<-as.data.frame(t(edata))
normlize<-normalize_input(as.matrix(data_lip))
reduce_UNTARRR<-Rtsne(normlize)$y
reduce_UNTARRR




data_untar<-as.data.frame(t(edata))
normlize<-normalize_input(as.matrix(data_untar))
reduce_UNTARRR<-Rtsne(normlize)$y
reduce_UNTARRR





####Incorporate additional PACE metadata 
add_data<-read.csv("Baker_PACE_data.csv",check.names=FALSE,na.strings=".")[,1:10]
add_data$Sex<-as.factor(add_data$Sex)
add_data$APOE<-as.factor(add_data$APOE)
final_join<-dplyr::inner_join(add_data,combined,by="StudyID")
write.csv(final_join,file="Wake_Lipidyzer.csv")










####Fit linear model for independent groups contrast
fit<-lmFit(edata,design,block=pheno_table$StudyID,correlation=corfit$consensus)
cm <- makeContrasts(
	`PACE2pre-PACEMCIpre` = PACE2pre-PACEMCIpre,
	`PACE2pre-PACEIGTpre` = PACE2pre-PACEIGTpre,
	`PACEMCIpre-PACEIGTpre` = PACEMCIpre-PACEIGTpre,
	`PACE2post-PACEMCIpost` = PACE2post-PACEMCIpost,
	`PACE2post-PACEIGTpost` = PACE2post-PACEIGTpost,
	`PACEMCIpost-PACEIGTpost` = PACEMCIpost-PACEIGTpost,
	levels=design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
T<-topTable(fit2,number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"FeatureID"

####Write DE table to file
Final_DE<-dplyr::inner_join(as.data.frame(featureData),T,by="FeatureID")%>%filter(Mode=="-")%>%select("MZ","RT","P.Value","F")
write.csv(Final_DE,file="Wake_Final_DE_Untargeted_BETWEEN_NEG.csv")
Final_DE<-dplyr::inner_join(as.data.frame(featureData),T,by="FeatureID")%>%filter(Mode=="+")%>%select("MZ","RT","P.Value","F")
write.csv(Final_DE,file="Wake_Final_DE_Untargeted_BETWEEN_POS.csv")






mummichog -f Wake_Final_DE_Untargeted_WITHIN_NEG_AEROBIC.txt -o Wake_WITHIN_NEG_AEROBIC -m negative -u 7 
mummichog -f Wake_Final_DE_Untargeted_WITHIN_POS.txt -o Wake_WITHIN_POS -m positive -u 7 -c .05
mummichog -f Wake_Final_DE_Untargeted_BETWEEN_NEG.txt -o Wake_BETWEEN_NEG -m negative -u 7
mummichog -f Wake_Final_DE_Untargeted_BETWEEN_POS.txt -o Wake_BETWEEN_POS -m positive -u 7