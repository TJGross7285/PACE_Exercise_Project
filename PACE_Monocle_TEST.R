library(monocle)




lipid<-read.csv("PACE_Lipidyzer_1_7.csv",check.names=FALSE)[,-c(1)]
energy<-read.csv("PACE_Polar_Metabolites_Energy_1_7.csv",check.names=FALSE)


energy_pheno<-energy[,1:5]
lipid_pheno<-lipid[,1:5]
lipid_meta<-as.data.frame(colnames(lipid[,-c(1:5)]))
colnames(lipid_meta)[1]<-"gene_short_name"
lipid_meta$gene_short_name<-as.character(lipid_meta$gene_short_name)
energy_meta<-as.data.frame(colnames(energy[,-c(1:5)]))
colnames(energy_meta)[1]<-"gene_short_name"

####Use lipid pheno as master pheno
all.equal(energy_pheno$StudyID,lipid_pheno$StudyID)

####Initialize metadata objects and total expression set
fd_lipid <- new("AnnotatedDataFrame", data = lipid_meta)
fd_energy<-new("AnnotatedDataFrame", data = energy_meta)
pd_lipid <- new("AnnotatedDataFrame", data = lipid_pheno)
pd_energy<-new("AnnotatedDataFrame", data = energy_pheno)


edata_lipid<-as.data.frame(t(lipid[,-c(1:5)]))
edata_energy<-as.data.frame(t(energy[,-c(1:5)]))


HSMM_lipid <- newCellDataSet(as.matrix(edata_lipid),
    phenoData = pd_lipid, featureData = fd_lipid, expressionFamily= uninormal())

HSMM_energy <- newCellDataSet(as.matrix(edata_energy),
    phenoData = pd, featureData = fd_energy, expressionFamily= uninormal())