library(limma)
library(dplyr)
library(impute)
library(sva)
library(GOplot)

####Read in and consolidate metadata from AKC lab 
A<-read.csv("A.csv",check.names=FALSE)
B<-read.csv("B.csv",check.names=FALSE)[,c(1:5)]

join<-inner_join(A,B,by="SampleID")
colnames(join)[15]<-"Key"
colnames(join)[14]<-"Order_Key_USE"
join$Order_Key_USE<-as.factor(join$Order_Key_USE)


write.csv(join,file="Wake_Pheno.csv")

####Import lipidyzer abundances
C<-read.csv("Lipid_Species_Concentration.csv",header=FALSE, na.strings=c(".","0"))
spec<-C[-c(1:5),-c(1:3)]
colnames(spec)<-unlist(C[4,4:dim(C)[2]])
D<-read.csv("Lipid_Class_Concentration.csv",header=FALSE, na.strings=c(".","0"))
class<-D[-c(1:5),-c(1:3)]
colnames(class)<-unlist(D[4,4:dim(D)[2]])
E<-read.csv("Fatty_Acid_Concentration.csv",header=FALSE, na.strings=c(".","0"))
FA<-E[-c(1:5),-c(1:3)]
colnames(FA)<-unlist(E[4,4:dim(E)[2]])
F<-read.csv("Total_Fatty_Acid.csv",header=FALSE, na.strings=c(".","0"))
total<-F[-c(1:5),-c(1:3)]
colnames(total)<-unlist(F[4,4:dim(F)[2]])

####Check that Georgetown IDs are same across tables 
identical(C[-c(1:5),1],D[-c(1:5),1],E[-c(1:5),1],F[-c(1:5),1])
Order_Key_USE<-as.factor(C[-c(1:5),1])

####Combine tables and join to metadata; threshold those with excessive NA 
bind<-cbind(Order_Key_USE,spec,class,FA,total)
combined<-inner_join(join,bind,by="Order_Key_USE")
na<-apply(combined,2,is.na)
na_index<-apply(na,2,sum)/dim(combined)[1]
combined<-combined[,na_index<.333]

####Incorporate additional PACE metadata 
add_data<-read.csv("Baker_PACE_data.csv",check.names=FALSE,na.strings=".")[,1:10]
add_data$Sex<-as.factor(add_data$Sex)
add_data$APOE<-as.factor(add_data$APOE)
final_join<-dplyr::inner_join(add_data,combined,by="StudyID")
levels(final_join$APOE)<-list("HasE4"=c("3_4","4_4","2_4"),"NoE4"=c("2_3","3_3"))
levels(final_join$Sex)<-list("F"=c("f","F"),"M"=c("M"))
write.csv(final_join,file="Wake_Lipidyzer.csv")

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

PACE2_exercise<-R%>%filter(Study=="PACE 2" & Treatment=="Aerobic")%>%select(-c(1:24))
PACE2_stretching<-R%>%filter(Study=="PACE 2" & Treatment=="Stretch")%>%select(-c(1:24))

###Check for the presence of metabolites with excessively missing/low variance abundances
gsg<-goodSamplesGenes(PACE2_exercise, verbose = 3)
gsg<-goodSamplesGenes(PACE2_stretching, verbose = 3)

### Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
###Call the network topology analysis function
sft <- pickSoftThreshold(PACE2_exercise, powerVector = powers, verbose = 5)








##1/27
red<-cbind(colnames(datExpr)[colors=="red"],rep("red",length(colnames(datExpr)[colors=="red"])))
green<-cbind(colnames(datExpr)[colors=="green"],rep("green",length(colnames(datExpr)[colors=="green"])))
black<-cbind(colnames(datExpr)[colors=="black"],rep("black",length(colnames(datExpr)[colors=="black"])))
grey<-cbind(colnames(datExpr)[colors=="grey"],rep("grey",length(colnames(datExpr)[colors=="grey"])))
turquoise<-cbind(colnames(datExpr)[colors=="turquoise"],rep("turquoise",length(colnames(datExpr)[colors=="turquoise"])))
yellow<-cbind(colnames(datExpr)[colors=="yellow"],rep("yellow",length(colnames(datExpr)[colors=="yellow"])))
brown<-cbind(colnames(datExpr)[colors=="brown"],rep("brown",length(colnames(datExpr)[colors=="brown"])))
blue<-cbind(colnames(datExpr)[colors=="blue"],rep("blue",length(colnames(datExpr)[colors=="blue"])))

total<-as.data.frame(rbind(red,green,black,grey,turquoise,yellow,brown,blue))
colnames(total)<-c("Feature","Eigenlipid")
total$Feature<-as.character(total$Feature)
final_frame_words<-total%>%unnest_tokens(word,Feature)%>%count(word,Eigenlipid, sort = TRUE)%>%bind_tf_idf(word, Eigenlipid, n)%>%cast_sparse(Eigenlipid,word, tf_idf)
final_frame_bigram<-total%>%unnest_tokens(word,Feature,token='ngrams',n=2)%>%count(word,Eigenlipid, sort = TRUE)%>%bind_tf_idf(word, Eigenlipid, n)%>%cast_sparse(Eigenlipid,word, tf_idf)
final_frame_trigram<-total%>%unnest_tokens(word,Feature,token='ngrams',n=3)%>%count(word,Eigenlipid, sort = TRUE)%>%bind_tf_idf(word, Eigenlipid, n)%>%cast_sparse(Eigenlipid,word, tf_idf)
shingle<-total%>%unnest_tokens(word,Feature,token="character_shingles")%>%count(word,Eigenlipid, sort = TRUE)%>%bind_tf_idf(word, Eigenlipid, n)%>%cast_sparse(Eigenlipid,word, tf_idf)

final_frame_total<-cbind(final_frame_words,final_frame_bigram,final_frame_trigram)%>%as.matrix%>%t()
write.csv(final_frame_total,file="Big_Data_1_27_MJ.csv")


write.csv(total,file="PACE_Lipidyzer_1_7.csv")

#####Sex+APOE+Bsln_BMI+Bsln_Cholesterol+Bsln_LDL+Bsln_HDL+Bsln_Trigly+Bsln_Fasting_Glu
design1<-model.matrix(~0+as.factor(Treat),data=total)
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

####Fit linear model for within subjects contrast
fit1<-lmFit(edata,design1,block=pheno_table$StudyID,correlation=corfit1$consensus)

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

fit4_F <- contrasts.fit(fit1, cm4)
fit4_F <- eBayes(fit4_F,trend=TRUE)

fit5_F<- contrasts.fit(fit1, cm5)
fit5_F<-eBayes(fit5_F,trend=TRUE)

fit6_F<-contrasts.fit(fit1, cm6)
fit6_F<-eBayes(fit6_F,trend=TRUE)

A<-topTableF(fit1_F,number=100000); A$FeatureID<-rownames(A)
A<-A%>%filter(P.Value<.05)%>%arrange(abs(desc(PACE2post.Aerobic.PACE2pre.Aerobic)))%>%dplyr::select(FeatureID)

B<-topTableF(fit2_F,number=100000); B$FeatureID<-rownames(B)
B<-B%>%filter(P.Value<.05)%>%arrange(abs(desc(PACEIGTpost.Aerobic.PACEIGTpre.Aerobic)))%>%dplyr::select(FeatureID)

C<-topTableF(fit3_F,number=100000); C$FeatureID<-rownames(C)
C<-C%>%filter(P.Value<.05)%>%arrange(abs(desc(PACEMCIpost.Aerobic.PACEMCIpre.Aerobic)))%>%dplyr::select(FeatureID)

D<-topTableF(fit4_F,number=100000); D$FeatureID<-rownames(D)
D<-D%>%filter(P.Value<.05)%>%arrange(abs(desc(PACE2post.Stretching.PACE2pre.Stretching)))%>%dplyr::select(FeatureID)

E<-topTableF(fit5_F,number=100000); E$FeatureID<-rownames(E)
E<-E%>%filter(P.Value<.05)%>%arrange(abs(desc(PACEIGTpost.Stretching.PACEIGTpre.Stretching)))%>%dplyr::select(FeatureID)

F<-topTableF(fit6_F,number=100000); F$FeatureID<-rownames(F)
F<-F%>%filter(P.Value<.05)%>%arrange(abs(desc(PACEMCIpost.Stretching.PACEMCIpre.Stretching)))%>%dplyr::select(FeatureID)

plotting<-cbind(as.data.frame(total),as.data.frame(t(edata)))
String<- numeric(ncol(plotting[,-c(1:5)])) 
for(i in 6:ncol(plotting)){
               String[i] <- bd.test(plotting[,i]~plotting$Treat)$p.value
}



String<- numeric(nrow(edata)) 
for(i in 1:nrow(edata)){
                String[i] <- bd.test(edata[i,]~pheno)$p.value
}

String<- numeric(nrow(t(MEs))) 
for(i in 1:nrow(t(MEs))){
                String[i] <- bd.test(t(MEs)[i,]~pheno)$p.value
}


pdf(file="Ring.pdf")
for(i in 1:length(A)){
	vioplot(edata[A[i],]~pheno)
}
dev.off()

C<-topTableF(fit3_F,number=100000)%>%filter(P.Value<.05)%>%arrange(abs(desc(PACEMCIpost.Aerobic.PACEMCIpre.Aerobic)))
D<-topTableF(fit4_F,number=100000)%>%filter(P.Value<.05)%>%arrange(abs(desc(PACE2post.Stretching.PACE2pre.Stretching)))
E<-topTableF(fit5_F,number=100000)%>%filter(P.Value<.05)%>%arrange(abs(desc(PACEIGTpost.Stretching.PACEIGTpre.Stretching)))
F<-topTableF(fit6_F,number=100000)%>%filter(P.Value<.05)%>%arrange(abs(desc(PACEMCIpost.Stretching.PACEMCIpre.Stretching)))


pdf(file="PACE_Volcano_12_23_LIPIDYZER.pdf")
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








plotting<-cbind(as.data.frame(total),as.data.frame(t(edata)))
pdf(file="PACE_ViolPlot_trimmmm.pdf") 




pdf(file="WithCE.pdf")
for (i in 15:ncol(plotting)){
	vioplot(plotting[,i] ~ plotting$Main, main=colnames(plotting[i]), horizontal=TRUE, ylab="Log2 Abundances (uM)")
}
dev.off()



String<- numeric(ncol(plotting[,-c(1:5)])) 
for(i in 6:ncol(plotting)){
               String[i] <- bd.test(plotting[,i]~plotting$Treat)$p.value
}










total<-cbind(Treat,R[,4:11])
mod<-model.matrix(~as.factor(Treat)+Sex+APOE+Bsln_BMI+Bsln_Cholesterol+Bsln_LDL+Bsln_HDL+Bsln_Trigly+Bsln_Fasting_Glu,data=total)
total<-cbind(as.data.frame(Treat),as.data.frame(R[,4:11]))
colnames(total)[c(1,10:18)]<-c("Main",
					"SV1",
					"SV2",
					"SV3",
					"SV4",
					"SV5",
					"SV6",
					"SV7",
					"SV8",
					"SV9")
####Set up accessory objects for DE incorporating pre-post timepoints
design1<-model.matrix(~0+Main+SV1+SV2+SV3+SV4+SV5+SV6+SV7+SV8+SV9,data=total)
colnames(design1)[1:4]<-c("Aerobic.pre",
						   "Aerobic.post",
						   "Stretching.pre",
						   "Stretching.post")
corfit1<-duplicateCorrelation(edata,design1,block=pheno_table$StudyID)



####Fit linear model for pairwise contrasts 
fit1<-lmFit(edata,design1,block=pheno_table$StudyID,correlation=corfit1$consensus)
cm <- makeContrasts(
	`Aerobic.post-Aerobic.pre` = Aerobic.post-Aerobic.pre,
	`Stretching.post-Stretching.pre` = Stretching.post-Stretching.pre,
	`Aerobic.pre-Stretching.pre` = Aerobic.pre-Stretching.pre,
	`Aerobic.post-Stretching.post` = Aerobic.post-Stretching.post,
	levels=design1)
fit1_F <- contrasts.fit(fit1, cm)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"

write.csv(T, file="PACE_DE_Collapsed_6_25.csv")


list1<-T%>%filter(P.Value<.05 & abs(Aerobic.post.Aerobic.pre)>.1)%>%select(Feature,Aerobic.post.Aerobic.pre)
list2<-T%>%filter(P.Value<.05 & abs(Stretching.post.Stretching.pre)>.1)%>%select(Feature,Stretching.post.Stretching.pre)
list3<-T%>%filter(P.Value<.05 & abs(Aerobic.pre.Stretching.pre)>.1)%>%select(Feature, Aerobic.pre.Stretching.pre)
list4<-T%>%filter(P.Value<.05 & abs(Aerobic.post.Stretching.post)>.1)%>%select(Feature,Aerobic.post.Stretching.post)



pdf(file="PACE_DE_Collapsed_Venn_6_25.pdf")

GOVenn(list1,list2,label=c("Aerobic.post-Aerobic.pre",
								"Stretching.post-Stretching.pre"),
								lfc.col=c("seagreen","gold","firebrick2"))

GOVenn(list3,list4,label=c("Aerobic.pre-Stretching.pre",
								"Aerobic.post-Stretching.post"),
								lfc.col=c("seagreen","gold","firebrick2"))

dev.off()
















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




























T<-topTable(fit1_F,coef=1,number=100000)
T_anno<-rep("PACEMCIpost.Aerobic-PACEMCIpre.Aerobic",dim(T)[1])
correctedFC<-predFCm(fit1_F,coef=1)
T<-cbind(T,correctedFC)

Y<-topTable(fit1_F,coef=2,number=100000)
Y_anno<-rep("PACEMCIpost.Stretching-PACEMCIpre.Stretching",dim(Y)[1])
correctedFC<-predFCm(fit1_F,coef=2)
Y<-cbind(Y,correctedFC)

U<-topTable(fit1_F,coef=3,number=100000)
U_anno<-rep("PACEMCIpre.Aerobic-PACEMCIpre.Stretching",dim(U)[1])
correctedFC<-predFCm(fit1_F,coef=3)
U<-cbind(U,correctedFC)

I<-topTable(fit1_F,coef=4,number=100000)
I_anno<-rep("PACEMCIpost.Aerobic-PACEMCIpost.Stretching",dim(I)[1])
correctedFC<-predFCm(fit1_F,coef=4)
I<-cbind(I,correctedFC)

####O<-topTable(fit1_F,coef=5,number=100000)
####O_anno<-rep("PACEIGTpost.Stretching-PACEIGTpre.Stretching",dim(O)[1])
####correctedFC<-predFCm(fit1_F,coef=5)
####O<-cbind(O,correctedFC)

####P<-topTable(fit1_F,coef=6,number=100000)
####P_anno<-rep("PACEMCIpost.Stretching-PACEMCIpre.Stretching",dim(P)[1])
####correctedFC<-predFCm(fit1_F,coef=6)
####P<-cbind(P,correctedFC)

study_vector<-c(T_anno,Y_anno,U_anno,I_anno)
feature_vector<-c(rownames(T),rownames(Y),rownames(U),rownames(I))
total_tables<-rbind(T,Y,U,I)
final_frame<-cbind(feature_vector,study_vector,total_tables)
colnames(final_frame)[1]<-"Feature"
colnames(final_frame)[2]<-"Contrast"

write.csv(final_frame, file="PACE_Lipidyzer_6_12_PRE_POST_Contrasts.csv")



list1<-final_frame%>%filter(Contrast=="PACEMCIpost.Aerobic-PACEMCIpre.Aerobic" & P.Value<.05 & abs(logFC)>.1)%>%select(Feature,logFC)
list2<-final_frame%>%filter(Contrast=="PACEMCIpost.Stretching-PACEMCIpre.Stretching" & P.Value<.05 & abs(logFC)>.1)%>%select(Feature,logFC)
list3<-final_frame%>%filter(Contrast=="PACEMCIpre.Aerobic-PACEMCIpre.Stretching" & P.Value<.05 & abs(logFC)>.1)%>%select(Feature,logFC)
list4<-final_frame%>%filter(Contrast=="PACEMCIpost.Aerobic-PACEMCIpost.Stretching" & P.Value<.05 & abs(logFC)>.1)%>%select(Feature,logFC)



boxplot(as.data.frame(t(edata))~)


levels(Treat)<-list("PACE2pre.Aerobic"=levels(Treat)[1],
						"PACE2pre.Stretch"=levels(Treat)[2],
						"PACE2post.Aerobic"=levels(Treat)[3],
						"PACE2post.Stretching"=levels(Treat)[4],
						"PACEIGTpre.Aerobic"=levels(Treat)[5],
						"PACEIGTpre.Stretching"=levels(Treat)[6],
						"PACEIGTpost.Aerobic"=levels(Treat)[7],
						"PACEIGTpost.Stretching"=levels(Treat)[8],
						"PACEMCIpre.Aerobic"=levels(Treat)[9],
						"PACEMCIpre.Stretching"=levels(Treat)[10],
						"PACEMCIpost.Aerobic"=levels(Treat)[11],
						"PACEMCIpost.Stretching"=levels(Treat)[12])
####Box plot final metabolite matrix to conduct a final check for metabolites with poor distribution





pdf(file="PACE_Boxplots_Lipidyzer_6_17.pdf") 
par(mar=c(15,4,4,2))

R<-read.csv("Wake_Lipidyzer.csv",check.names=FALSE)[-c(135,136),]
pheno_table<-R[,c(2,3,14,18)]
colnames(pheno_table)<-c("StudyID","Treatment","Study","Main")

edata<-log2(t(R[,-c(1:25)]))

Treat<-factor(paste(pheno_table$Study,pheno_table$Main,sep="."))
Treat<-factor(paste(Treat,pheno_table$Treatment,sep="."))
levels(Treat)<-list("PACE2pre.Aerobic"=levels(Treat)[1],
						"PACE2pre.Stretch"=levels(Treat)[2],
						"PACE2post.Aerobic"=levels(Treat)[3],
						"PACE2post.Stretching"=levels(Treat)[4],
						"PACEIGTpre.Aerobic"=levels(Treat)[5],
						"PACEIGTpre.Stretching"=levels(Treat)[6],
						"PACEIGTpost.Aerobic"=levels(Treat)[7],
						"PACEIGTpost.Stretching"=levels(Treat)[8],
						"PACEMCIpre.Aerobic"=levels(Treat)[9],
						"PACEMCIpre.Stretching"=levels(Treat)[10],
						"PACEMCIpost.Aerobic"=levels(Treat)[11],
						"PACEMCIpost.Stretching"=levels(Treat)[12])

D<-cbind(Treat,as.data.frame(t(edata)))
D<-D%>%filter(Treat=="PACEMCIpre.Aerobic"|Treat=="PACEMCIpre.Stretching"|Treat=="PACEMCIpost.Aerobic"|Treat=="PACEMCIpost.Stretching")
Treat<-D$Treat
edata<-D[,-c(1)]

for (i in 1:ncol(edata)){
	boxplot(edata[,i] ~ factor(Treat), main=colnames(edata)[i],ylab= "Log2(uM)",xlab=NULL, las=2)
}
dev.off()


vec<-vector()
for (i in 1:ncol(edata)){
	inter<-ad.test(edata[,i]~factor(Treat))
	vec[i]<-inter$ad[1,3]
}
names(vec)<-colnames(edata)







































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






####Define custom sum
Treat<-factor(paste(pheno_table$Study,pheno_table$Treatment,sep="."))
total<-cbind(Treat,R[,4:11])
design1<-model.matrix(~0+Treat+Sex+APOE+Bsln_BMI+Bsln_Cholesterol+Bsln_LDL+Bsln_HDL+Bsln_Trigly+Bsln_Fasting_Glu,data=total)
colnames(design1)[1:6]<-c("PACE2.Aerobic",
						"PACE2.Stretching",
						"PACEIGT.Aerobic",
						"PACEIGT.Stretching",
						"PACEMCI.Aerobic",
						"PACEMCI.Stretching")
corfit1<-duplicateCorrelation(edata,design1,block=pheno_table$StudyID)

frame<-as.data.frame(t(edata))
select<-R%>%filter(Treat=="")

####Fit linear model for within subjects contrast
fit1<-lmFit(edata,design1,block=pheno_table$StudyID,correlation=corfit1$consensus)
cm <- makeContrasts(
	`PACE2post.Aerobic-PACE2pre.Aerobic` = PACE2post.Aerobic-PACE2pre.Aerobic,
	`PACEIGTpost.Aerobic-PACEIGTpre.Aerobic` = PACEIGTpost.Aerobic-PACEIGTpre.Aerobic,
	`PACEMCIpost.Aerobic-PACEMCIpre.Aerobic` = PACEMCIpost.Aerobic-PACEMCIpre.Aerobic,
	`PACE2post.Stretching-PACE2pre.Stretching` = PACE2post.Stretching-PACE2pre.Stretching,
	`PACEIGTpost.Stretching-PACEIGTpre.Stretching` = PACEIGTpost.Stretching-PACEIGTpre.Stretching,
	`PACEMCIpost.Stretching-PACEMCIpre.Stretching` = PACEMCIpost.Stretching-PACEMCIpre.Stretching,
	`pace
	levels=design1)
fit1_F <- contrasts.fit(fit1, cm)
fit1_F <- eBayes(fit1_F,trend=TRUE,robust=TRUE)

























write.csv(T,file="Wake_Final_DE_Lipidyzer.csv")

####Calculate correlations of metabolites to other pheno info 
library(gplots)
R<-read.csv("Wake_Lipidyzer.csv",check.names=FALSE)
pheno_corr<-apply(R[,c(6:11)],2,as.numeric)
abunds<-apply(R[,-c(1:25)],2,as.numeric)
corr<-cor(abunds,pheno_corr,method="spearman",use="complete.obs")
pdf(file="Lipidyzer_Pheno_Corr.pdf")
heatmap.2(corr,trace="none")
dev.off()






























toCluster<-T%>%filter(P.Value<.05)
toCluster<-t(T[,2:7])
####Write DE table to file
write.csv(T,file="Wake_Final_DE_Lipidyzer.csv")






fit<-lmFit(edata,design,block=pheno_table$StudyID,correlation=corfit$consensus)
cm <- makeContrasts(
	`PACE2post.Stretching-PACE2pre.Stretching` = PACE2post.Stretching-PACE2pre.Stretching,
	`PACEIGTpost.Stretching-PACEIGTpre.Stretching` = PACEIGTpost.Aerobic-PACEIGTpre.Aerobic,
	`PACEMCIpost.Stretching-PACEMCIpre.Stretching` = PACEMCIpost.Stretching-PACEMCIpre.Stretching,
	`PACE2pre.Stretching-PACEMCIpre.Stretching`= PACE2pre.Aerobic-PACEMCIpre.Aerobic,
	`PACE2pre.Stretching-PACEIGTpre.Stretching`= PACE2pre.Aerobic-PACEIGTpre.Aerobic,
	`PACEMCIpre.Stretching-PACEIGTpre.Stretching`= PACEMCIpre.Aerobic-PACEIGTpre.Aerobic,
	`PACE2post.Stretching-PACEMCIpost.Stretching`= PACE2pre.Aerobic-PACEMCIpre.Aerobic,
	`PACE2post.Stretching-PACEIGTpost.Stretching`= PACE2pre.Aerobic-PACEIGTpre.Aerobic,
	`PACEMCIpost.Stretching-PACEIGTpost.Stretching`= PACEMCIpre.Aerobic-PACEIGTpre.Aerobic,
	levels=design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2,robust=TRUE)
T<-topTable(fit2,number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"FeatureID"

write.csv(T,file="Wake_Final_DE_Lipidyzer_CONTROL.csv")






####Box plot final metabolite matrix to conduct a final check for metabolites with poor distribution
pdf(file="Lipidyzer_WAKE_Boxplot.pdf") 
plot<-plot[pheno_table$Treatment=="Stretch",]
pheno<-pheno_table$Main[pheno_table$Treatment=="Stretch"]
for (i in 1:ncol(plot)){
	boxplot(plot[, i] ~ pheno, main=colnames(plot[i]))
}
dev.off()




R<-cbind(as.data.frame(pheno_table),as.data.frame(t(edata)))

A<-R%>%filter(Study=="PACE 2" & Main=="1" & Treatment=="Aerobic")%>%select("FFA(18:2)")
summary(A$`FFA(18:2)`)


R<-cbind(as.data.frame(pheno_table),as.data.frame(t(edata)))

B<-R%>%filter(Study=="PACE 2" & Main=="2" & Treatment=="Aerobic")%>%select("FFA(18:2)")
summary(B$`FFA(18:2)`)


E<-model.matrix()



















####Fit linear model
fit<-lmFit(edata,design,block=pheno_table$StudyID,correlation=corfit$consensus)
cm <- makeContrasts(
	`PACE2pre-PACE2post` = PACE2pre-PACE2post,
	`PACEIGTpre-PACEIGTpost` = PACEIGTpre-PACEIGTpost,
	`PACEMCIpre-PACEMCIpost` = PACEMCIpre-PACEMCIpost,
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
Final_DE<-dplyr::inner_join(as.data.frame(featureData),T,by="FeatureID")
write.csv(Final_DE,file="Wake_Final_DE_Untargeted.csv")


























combined<-cbind(as.data.frame(R[,c(2:16)]),as.data.frame(t(edata)))

Pace2_abunds<-t(combined%>%filter(Study=="PACE 2")%>%select(-c(1:15)))
Pace2_pheno<-combined%>%filter(Study=="PACE 2")%>%select(`Visit #`)
Pace2_model<-model.matrix(~as.factor(`Visit #`),data=Pace2_pheno)

PaceIGT_abunds<-t(combined%>%filter(Study=="PACE-IGT")%>%select(-c(1:15)))
PaceIGT_pheno<-combined%>%filter(Study=="PACE-IGT")%>%select(`Visit #`)
PaceIGT_model<-model.matrix(~as.factor(`Visit #`),data=PaceIGT_pheno)

PaceMCI_abunds<-t(combined%>%filter(Study=="PACE-MCI")%>%select(-c(1:15)))
PaceMCI_pheno<-combined%>%filter(Study=="PACE-MCI")%>%select(`Visit #`)
PaceMCI_model<-model.matrix(~as.factor(`Visit #`),data=PaceMCI_pheno)

library(limma)
fit_pace2<-lmFit(Pace2_abunds,Pace2_model)
fit_pace2<-eBayes(fit_pace2)
pace2<-topTable(fit_pace2, adjust="BH",number=10000)
write.csv(pace2,file="Wake_DE_PACE_2.csv")

fit_paceIGT<-lmFit(PaceIGT_abunds,PaceIGT_model)
fit_paceIGT<-eBayes(fit_paceIGT)
paceIGT<-topTable(fit_paceIGT, adjust="BH",number=10000)
write.csv(paceIGT,file="Wake_DE_PACE_IGT.csv")

fit_paceMCI<-lmFit(PaceMCI_abunds,PaceMCI_model)
fit_paceMCI<-eBayes(fit_paceMCI)
paceMCI<-topTable(fit_paceMCI, adjust="BH",number=10000)
write.csv(paceMCI,file="Wake_DE_PACE_MCI.csv")































####Calculate SVs 
library(sva)
Dementia<-R$Study
levels(Dementia)<-list("Dementia"=c("PACE-MCI","PACE 2"),"Non-Demented"=c("PACE-IGT"))
Diabetes<-R$Study
levels(Diabetes)<-list("Diabetes"=c("PACE-IGT","PACE 2"),"Non-Diabetes"=c("PACE-MCI"))
design_table<-cbind(as.data.frame(as.factor(R$`Visit #`)),as.data.frame(as.factor(R$`Study`)),as.data.frame(as.factor(Dementia)),as.data.frame(as.factor(Diabetes)))

mod<-model.matrix(~Dementia+Diabetes,data=design_table)
mod0<-model.matrix(~1,data=as.factor(R$`Visit #`))
n.sv<-num.sv(edata,mod,seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)


total_table<-cbind(design_table,as.data.frame(svobj$sv))
colnames(total_table)<-c("visit","study","dementia","diabetes","SV1","SV2","SV3","SV4","SV5","SV6","SV7","SV8")
model<-model.matrix(~dementia:diabetes:visit+SV1+SV2+SV3+SV4+SV5+SV6+SV6+SV7+SV8, data=total_table)
model<-model[,-c(13,17)]
model<-model[,-c(10:12)]

library(limma)
fit<-lmFit(edata,model)
cm <- makeContrasts(
+ PACE_2 = `PACE 2.2`-`PACE 2.1`,
+ PACE_IGT = `PACE-IGT.2`-`PACE-IGT.1`,
+ PACE_MCI = `PACE-MCI.2`-`PACE-MCI.1`,
+ levels=design)
fit<-eBayes(fit)
T<-topTable(fit, adjust="BH",number=10000)
write.csv(T,file="Wake_DE.csv")