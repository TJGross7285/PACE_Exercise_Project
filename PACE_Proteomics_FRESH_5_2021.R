library(dplyr)

setwd("/Volumes/NO\ NAME/Active_Projects_4_9_2021/PACE\ Proteomics\ ")
A<-read.csv("Clean_PACE_Prot_9_24.csv",check.names=FALSE)[,-c(20:21)]%>%filter(Study=="PACE-MCI")%>%select(-c("PPY.1","PPY.2"))

pheno<-A[,c(3:11)]

abunds_1<-A[,25:43]
abunds_2<-A[,44:62]

diff<-(abunds_2-abunds_1)/(abunds_2+abunds_1)

##########INTERESTING FEATURES: @@@3(treat,HDL trend), 5(Bsln_Trigly), 6(Bsln_HDL), 9(Sex trend), @@@@@11(Treat,tot Cholest) IL-6, 12 (Treat trend) IL-7,
########## 13(APOE trend), 15 (Sex), @@@@16(Treat,LDL) TARC, @@@@18(Treat,Sex) TNF-A
anova(lm(diff[,1]~as.factor(pheno$Treatment)
	             +as.factor(pheno$Sex)
	             +as.factor(pheno$APOE)
	             +as.numeric(pheno$Bsln_BMI)
	             +as.numeric(pheno$Bsln_Cholesterol)
	             +as.numeric(pheno$Bsln_LDL)
	             +as.numeric(pheno$Bsln_HDL)
	             +as.numeric(pheno$Bsln_Trigly)
	             +as.numeric(pheno$Bsln_Fasting_Glu)))

vioplot(log2(dv[,3])~as.factor(pheno$Treatment))
vioplot(log2(dv[,11])~as.factor(pheno$Treatment))
vioplot(log2(dv[,16])~as.factor(pheno$Treatment))
vioplot(log2(dv[,18])~as.factor(pheno$Treatment))


p <- ggplot(ToothGrowth, aes(x=Treatment, y=len)) + 
  geom_violin()