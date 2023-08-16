library(dplyr)
library(impute)

A<-read.csv("Complete_Wake_Pheno.csv")[,-c(1)]
B<-read.csv("20200624_pace_lipidomics data.NIST.csv",check.names=FALSE)[,-c(16)]
B<-B%>%mutate(CV = gsub("%?","",B$`Coefficient of variation`))%>%filter(CV<=40)

E<-read.csv("20200624_pace_lipidomics data.Abunds.csv",check.names=FALSE)
Eprime<-E%>%t()
colnames(Eprime)<-Eprime[1,]
SampleID<-colnames(E)
bound<-cbind(as.data.frame(SampleID),as.data.frame(Eprime))
abunds<-bound[2:dim(bound)[1],-c(1)]
rownames(abunds)<-NULL

F<-read.csv("20200624_pace_lipidomics data.PoolQC.csv",check.names=FALSE)[,-c(25)]
F<-F%>%mutate(CV = gsub("%?","",F$`Coefficient of variation`))%>%filter(CV<=40)
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
lipid_bounds<-reshape(lipid_bound,timevar="Visit..",idvar="StudyID",direction="wide")[,-c(306:329)]







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
metab_bounds<-reshape(metab_bound,timevar="Visit..",idvar="StudyID",direction="wide")[,-c(121:144)]



all.equal(metab_bound$StudyID,lipid_bound$StudyID)



A<-read.csv("Complete_Wake_Pheno.csv")[,-c(1)]
W<-read.csv("PACE_Proteomics_7_2020.csv",check.names=FALSE)
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
prot_bounds<-reshape(prot_bound,timevar="Visit..",idvar="StudyID.x",direction="wide")[,-c(51:79)]


write.csv(lipid_bounds,file="Clean_PACE_Lipid_9_30_CV40.csv")
write.csv(metab_bounds,file="Clean_PACE_Metab_9_30_CV40.csv")
write.csv(prot_bounds,file="Clean_PACE_Prot_9_30_CV40.csv")

###########################################
###########################################
###########################################
###########################################
###########################################
###########################################
###########################################

library(dplyr)
library(impute)

A<-read.csv("Complete_Wake_Pheno.csv")[,-c(1)]
B<-read.csv("20200624_pace_lipidomics data.NIST.csv",check.names=FALSE)[,-c(16)]
B<-B%>%mutate(CV = gsub("%?","",B$`Coefficient of variation`))%>%filter(CV<=20)

E<-read.csv("20200624_pace_lipidomics data.Abunds.csv",check.names=FALSE)
Eprime<-E%>%t()
colnames(Eprime)<-Eprime[1,]
SampleID<-colnames(E)
bound<-cbind(as.data.frame(SampleID),as.data.frame(Eprime))
abunds<-bound[2:dim(bound)[1],-c(1)]
rownames(abunds)<-NULL

F<-read.csv("20200624_pace_lipidomics data.PoolQC.csv",check.names=FALSE)[,-c(25)]
F<-F%>%mutate(CV = gsub("%?","",F$`Coefficient of variation`))%>%filter(CV<=20)
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
lipid_bounds<-reshape(lipid_bound,timevar="Visit..",idvar="StudyID",direction="wide")[,-c(178:201)]







G<-read.csv("20200625_pace_metabolomics.NIST.csv",check.names=FALSE)[,-c(3,18)]
G<-G%>%mutate(CV = gsub("%?","",G$`Coefficient of variation`))%>%filter(CV<=20)

C<-read.csv("20200625_pace_metabolomics.Abunds.csv",check.names=FALSE)
Cprime<-C%>%t()
colnames(Cprime)<-Cprime[1,]
SampleID<-colnames(C)
bound<-cbind(as.data.frame(SampleID),as.data.frame(Cprime))
abunds<-bound[4:dim(bound)[1],-c(1)]
rownames(abunds)<-NULL

D<-read.csv("20200625_pace_metabolomics.PoolQC.csv",check.names=FALSE)[,-c(3,18)]
D<-D%>%mutate(CV = gsub("%?","",D$`Coefficient of variation`))%>%filter(CV<=20)
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
metab_bounds<-reshape(metab_bound,timevar="Visit..",idvar="StudyID",direction="wide")[,-c(72:95)]



all.equal(metab_bound$StudyID,lipid_bound$StudyID)



A<-read.csv("Complete_Wake_Pheno.csv")[,-c(1)]
W<-read.csv("PACE_Proteomics_7_2020.csv",check.names=FALSE)
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
prot_bounds<-reshape(prot_bound,timevar="Visit..",idvar="StudyID.x",direction="wide")[,-c(51:79)]


write.csv(lipid_bounds,file="Clean_PACE_Lipid_9_30_CV20.csv")
write.csv(metab_bounds,file="Clean_PACE_Metab_9_30_CV20.csv")


