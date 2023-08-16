A<-read.csv("20190830_Wake Forest Metabolomics Data.csv",check.names=FALSE,na.strings=".")
B<-read.csv("Complete_Wake_Pheno.csv",check.names=FALSE)[,-c(1)]
U<-dplyr::inner_join(B,A,by="EnergyID")
index<-duplicated(U$EnergyID)
write.csv(U[index==FALSE,],file="Energy_PACE_Fixed.csv")













R<-read.csv("Energy_PACE_Fixed.csv",check.names=FALSE)

pheno_table<-R[,c(2:29)]
pheno_index<-apply(apply(pheno_table,2,is.na),1,sum)
pheno_table<-pheno_table[pheno_index<1,]


edata<-R[pheno_index<1,-c(1:29)]
edata_index<-apply(apply(edata,2,is.na),2,sum)/dim(edata)[1]
edata<-log2(t(edata[,edata_index<.333]))
edata<-impute.knn(edata,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data

Treat<-factor(paste(pheno_table$Treatment,pheno_table$`Visit #`,sep="."))
Treat<-factor(paste(pheno_table$Study,Treat,sep="."))
Treat<-gsub("-","",Treat)
Treat<-gsub(" ","",Treat)
total<-cbind(Treat,pheno_table)
#########+Sex+APOE+Bsln_BMI
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

A<-topTableF(fit1_F,number=100000)
A<-cbind(rownames(A),A)
colnames(A)[1]<-"Feature"
A%>%filter(P.Value<.05)%>%arrange(abs(PACE2post.Aerobic.PACE2pre.Aerobic))

B<-topTableF(fit2_F,number=100000)
B<-cbind(rownames(B),B)
colnames(B)[1]<-"Feature"
B%>%filter(P.Value<.05)%>%arrange(desc(abs(PACEIGTpost.Aerobic.PACEIGTpre.Aerobic)))

C<-topTableF(fit3_F,number=100000)
C<-cbind(rownames(C),C)
colnames(C)[1]<-"Feature"

D<-topTableF(fit4_F,number=100000)
D<-cbind(rownames(D),D)
colnames(D)[1]<-"Feature"

E<-topTableF(fit5_F,number=100000)
E<-cbind(rownames(E),E)
colnames(E)[1]<-"Feature"

F<-topTableF(fit6_F,number=100000)
F<-cbind(rownames(F),F)
colnames(F)[1]<-"Feature"


pdf(file="PACE_Volcano_12_27_ENERGY.pdf")
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



















svobj<-sva(edata,mod,mod0,n.sv=n.sv)$sv
total<-cbind(total,as.data.frame(svobj))
colnames(total)[c(1,30:36)]<-c("Main",
					"SV1",
					"SV2",
					"SV3",
					"SV4",
					"SV5",
					"SV6",
					"SV7")

pheno_table<-cbind()

####Set up accessory objects for DE incorporating pre-post timepoints
design1<-model.matrix(~0+Main+SV1+SV2+SV3+SV4+SV5+SV6+SV7,data=total)
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

write.csv(T, file="PACE_DE_EnergyMetabolism_9_4.csv")
