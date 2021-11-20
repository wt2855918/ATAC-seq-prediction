[1]GM12878
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/ATAC_seq/Batch_1/ENCFF172DEA.bed"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/DNase_seq/ENCFF235KUD.bed"
myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/FAIRE-seq/ENCFF001UYE.bed"

res1 = read.table(myinf1,sep="\t",quote=NULL)
res2 = read.table(myinf2,sep="\t",quote=NULL)
res3 = read.table(myinf3,sep="\t",quote=NULL)

xx1 = res1$V3 - res1$V2
xx2 = res2$V3 - res2$V2
xx3 = res3$V3 - res3$V2

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_8/GM12878_ATAC_distribution.pdf"
pdf(myoutf,width=5,height=5)
plot(density(xx1))
dev.off()

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_8/GM12878_DNase_distribution.pdf"
pdf(myoutf,width=5,height=5)
plot(density(xx2))
dev.off()

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_8/GM12878_FAIRE_distribution.pdf"
pdf(myoutf,width=5,height=5)
plot(density(xx3))
dev.off()

[2]HepG2
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/ATAC_seq/ENCFF356TXH.bed"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/DNase_seq/ENCFF422EDI.bed"
myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/FAIRE-seq/ENCFF001UYN.bed"

res1 = read.table(myinf1,sep="\t",quote=NULL)
res2 = read.table(myinf2,sep="\t",quote=NULL)
res3 = read.table(myinf3,sep="\t",quote=NULL)

xx1 = res1$V3 - res1$V2
xx2 = res2$V3 - res2$V2
xx3 = res3$V3 - res3$V2

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_8/HepG2_ATAC_distribution.pdf"
pdf(myoutf,width=5,height=5)
plot(density(xx1))
dev.off()

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_8/HepG2_DNase_distribution.pdf"
pdf(myoutf,width=5,height=5)
plot(density(xx2))
dev.off()

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_8/HepG2_FAIRE_distribution.pdf"
pdf(myoutf,width=5,height=5)
plot(density(xx3))
dev.off()