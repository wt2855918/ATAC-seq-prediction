[1]GM12878
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/Histone_file_annotation_bigwig.txt"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/TF_file_annotation_bigWig.txt"
myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/Histone_file_annotation_bigwig.txt"
myinf4 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/TF_file_annotation_bigWig.txt"

res1 = read.table(myinf1,sep="\t",quote=NULL)
res2 = read.table(myinf2,sep="\t",quote=NULL)
res3 = read.table(myinf3,sep="\t",quote=NULL)
res4 = read.table(myinf4,sep="\t",quote=NULL)

info1 = rbind(res1,res2)
info2 = rbind(res1,res2)

info1$Cell = rep("GM12878",nrow(info1))
info2$Cell = rep("HepG2",nrow(info2))

info = rbind(info1,info2)

myoutf = "/ihome/yanding/ATAC_Suppl_Table_1.xls"
write.table(info,myoutf,sep="\t",quote=F)