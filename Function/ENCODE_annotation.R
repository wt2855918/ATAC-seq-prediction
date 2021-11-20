[1]GM12878 Annotation Histone
rm(list=ls())
myinf= "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/Histone_ChIP_seq_bigWig/metadata.tsv"
info = read.table(myinf,sep="\t",quote=NULL,header=T)
info = info[,c("File.accession","Experiment.target")]
info$Experiment.target = gsub("-human","",info$Experiment.target)
row.names(info) = info$File.accession

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/Histone_file_annotation_bigwig.txt"
write.table(info,myoutf,sep="\t",quote=F)

[2]GM12878 Annotation TF
rm(list=ls())
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/metadata.tsv"
info = read.table(myinf,sep="\t",quote=NULL,header=T)
info = info[,c("File.accession","Experiment.target")]
info$Experiment.target = gsub("-human","",info$Experiment.target)
row.names(info) = info$File.accession

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/TF_file_annotation_bigWig.txt"
write.table(info,myoutf,sep="\t",quote=F)

[3]HepG2 Annotation Histone
rm(list=ls())
myinf= "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/Histone_ChIP_seq_bigWig/metadata.tsv"
info = read.table(myinf,sep="\t",quote=NULL,header=T)
info = info[,c("File.accession","Experiment.target")]
info$Experiment.target = gsub("-human","",info$Experiment.target)
row.names(info) = info$File.accession

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/Histone_file_annotation_bigWig.txt"
write.table(info,myoutf,sep="\t",quote=F)

[4]HepG2 Annotation TF
rm(list=ls())
myinf= "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/metadata.tsv"
info = read.table(myinf,sep="\t",quote=NULL,header=T)
info = info[,c("File.accession","Experiment.target")]
info$Experiment.target = gsub("-human","",info$Experiment.target)
row.names(info) = info$File.accession

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/TF_file_annotation_bigWig.txt"
write.table(info,myoutf,sep="\t",quote=F)



[4]mESC Annotation Histone
rm(list=ls())
myinf= "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/E14TG2a/Histone_ChIP_seq_bigWig/metadata.tsv"
info = read.table(myinf,sep="\t",quote=NULL,header=T)
info = info[,c("File.accession","Experiment.target")]
info$Experiment.target = gsub("-mouse","",info$Experiment.target)
row.names(info) = info$File.accession

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/E14TG2a/Histone_file_annotation_bigWig.txt"

write.table(info,myoutf,sep="\t",quote=F)
