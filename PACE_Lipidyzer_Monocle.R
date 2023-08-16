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

lipid_Stretch<-cbind(as.data.frame(total),as.data.frame(t(edata)))%>%filter(Study=="PACE 2"&Treatment=="Stretch")
lipid_Aerobic<-cbind(as.data.frame(total),as.data.frame(t(edata)))%>%filter(Study=="PACE 2"&Treatment=="Aerobic")

stretch_abunds<-t(lipid_Stretch[,-c(1:5)])
aerobic_abunds<-t(lipid_Aerobic[,-c(1:5)])

stretch_pheno<-lipid_Stretch[,1:5]
aerobic_pheno<-lipid_Aerobic[,1:5]

stretch_features<-as.data.frame(rownames(stretch_abunds))
colnames(stretch_features)[1]<-"gene_short_name"
rownames(stretch_features)<-rownames(stretch_abunds)

aerobic_features<-as.data.frame(rownames(aerobic_abunds))
colnames(aerobic_features)[1]<-"gene_short_name"
rownames(aerobic_features)<-rownames(aerobic_abunds)

####Initialize metadata objects and total expression set
fd_lipid_aerobic<-new("AnnotatedDataFrame", data = aerobic_features)
pd_lipid_aerobic<-new("AnnotatedDataFrame", data = aerobic_pheno)
fd_lipid_stretch<-new("AnnotatedDataFrame", data = stretch_features)
pd_lipid_stretch<-new("AnnotatedDataFrame", data = stretch_pheno)


HSMM_lipid_aerobic <- newCellDataSet(as.matrix(aerobic_abunds),
    phenoData = pd_lipid_aerobic, featureData = fd_lipid_aerobic, expressionFamily= uninormal())

HSMM_lipid_stretch <- newCellDataSet(as.matrix(stretch_abunds),
    phenoData = pd_lipid_stretch, featureData = fd_lipid_stretch, expressionFamily= uninormal())

####Carry out dimension reduction/sample ordering using SVs as residual medel formula   
HSMM_myo_aerobic_lipid<- reduceDimension(HSMM_lipid_aerobic, reduction_method = 'DDRTree',norm_method="none", pseudo_expr=0)
HSMM_myo1_aerobic_lipid <- orderCells(HSMM_myo_aerobic_lipid)
plot_cell_trajectory(HSMM_myo1_aerobic_lipid,color_by = "Treat")


####Carry out dimension reduction/sample ordering using SVs as residual medel formula   
HSMM_myo_stretch_lipid<- reduceDimension(HSMM_lipid_stretch, reduction_method = 'DDRTree',norm_method="none", pseudo_expr=0)
HSMM_myo1_stretch_lipid<- orderCells(HSMM_myo_stretch_lipid)
plot_cell_trajectory(HSMM_myo1_stretch_lipid,color_by = "Treat")




######################################
######################################
######################################

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

energy_Stretch<-cbind(as.data.frame(total),as.data.frame(t(edata)))%>%filter(Study=="PACE 2"&Treatment=="Stretch")
energy_Aerobic<-cbind(as.data.frame(total),as.data.frame(t(edata)))%>%filter(Study=="PACE 2"&Treatment=="Aerobic")

stretch_abunds<-t(energy_Stretch[,-c(1:29)])
aerobic_abunds<-t(energy_Aerobic[,-c(1:29)])

stretch_pheno<-energy_Stretch[,1:29]
aerobic_pheno<-energy_Aerobic[,1:29]

stretch_features<-as.data.frame(rownames(stretch_abunds))
colnames(stretch_features)[1]<-"gene_short_name"
rownames(stretch_features)<-rownames(stretch_abunds)

aerobic_features<-as.data.frame(rownames(aerobic_abunds))
colnames(aerobic_features)[1]<-"gene_short_name"
rownames(aerobic_features)<-rownames(aerobic_abunds)

####Initialize metadata objects and total expression set
fd_lipid_aerobic<-new("AnnotatedDataFrame", data = aerobic_features)
pd_lipid_aerobic<-new("AnnotatedDataFrame", data = aerobic_pheno)
fd_lipid_stretch<-new("AnnotatedDataFrame", data = stretch_features)
pd_lipid_stretch<-new("AnnotatedDataFrame", data = stretch_pheno)


HSMM_energy_aerobic <- newCellDataSet(as.matrix(aerobic_abunds),
    phenoData = pd_lipid_aerobic, featureData = fd_lipid_aerobic, expressionFamily= uninormal())

HSMM_energy_stretch <- newCellDataSet(as.matrix(stretch_abunds),
    phenoData = pd_lipid_stretch, featureData = fd_lipid_stretch, expressionFamily= uninormal())

####Carry out dimension reduction/sample ordering using SVs as residual medel formula   
HSMM_myo_aerobic_energy<- reduceDimension(HSMM_energy_aerobic, reduction_method = 'DDRTree',norm_method="none", pseudo_expr=0)
HSMM_myo1_aerobic_energy<- orderCells(HSMM_myo_aerobic_energy)
plot_cell_trajectory(HSMM_myo1_aerobic_energy,color_by = "Treat")


####Carry out dimension reduction/sample ordering using SVs as residual medel formula   
HSMM_myo_stretch_energy<- reduceDimension(HSMM_energy_stretch, reduction_method = 'DDRTree',norm_method="none", pseudo_expr=0)
HSMM_myo1_stretch_energy<- orderCells(HSMM_myo_stretch_energy)
plot_cell_trajectory(HSMM_myo1_stretch_energy,color_by = "Treat")






pdf(file="PACE_Monocle_Trajectories_1_13.pdf")
plot_cell_trajectory(HSMM_myo1_aerobic_lipid,color_by = "Treat")
plot_cell_trajectory(HSMM_myo1_stretch_lipid,color_by = "Treat")
plot_cell_trajectory(HSMM_myo1_aerobic_energy,color_by = "Treat")
plot_cell_trajectory(HSMM_myo1_stretch_energy,color_by = "Treat")
dev.off()







