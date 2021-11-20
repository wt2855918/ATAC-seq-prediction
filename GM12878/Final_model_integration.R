library(randomForest)
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Multi_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Uni_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/moveme.R")

#########################################################
#load TF chip-seq profile
#########################################################
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_counts_bins_logFC/"
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/TF_file_annotation_bigWig.txt"

files = list.files(mydir)
name = gsub(".Rda","",files)
name = gsub("_ENCFF172DEA_positive","",name)
name = gsub("_ENCFF172DEA_negative","",name)
name = unique(name)

info = read.table(myinf1,sep="\t",quote=NULL,stringsAsFactors=F)
info = info[intersect(name,row.names(info)),]
colnames(info)= c("target","file_accession")

library(randomForest)
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Multi_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Uni_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/moveme.R")

AUC_score = list()

#load the positive and negative data
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_positive.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_negative.Rda")

#load the tag
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_neg.Rda")

pos = validation_pos
neg = validation_neg

#tag_pos = sample(seq(1,895736),5000)
#tag_neg = sample(seq(1,895736),5000)

rm(positive)
rm(negative)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_counts_bins_logFC/"
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/TF_file_annotation_bigWig.txt"

files = list.files(mydir)
name = gsub(".Rda","",files)
name = gsub("_ENCFF172DEA_positive","",name)
name = gsub("_ENCFF172DEA_negative","",name)
name = unique(name)

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/TF_file_annotation_bigWig.txt"
info = read.table(myinf1,sep="\t",quote=NULL,stringsAsFactors=F)
info = info[intersect(name,row.names(info)),]
colnames(info)= c("target","file_accession")

res = matrix(0, 4000, nrow(info))
colnames(res) = row.names(info)
res = as.data.frame(res)

for(i in 1 : nrow(info))
{
	cat("\r",i)
	
	tmpinf1 = paste0(mydir,row.names(info)[i],"_ENCFF172DEA_negative.Rda")
	load(tmpinf1)
	negative = histone
	negative = negative[row.names(neg),]
	
	tmpinf2 = paste0(mydir,row.names(info)[i],"_ENCFF172DEA_positive.Rda")
	load(tmpinf2)
	positive = histone
	positive = positive[row.names(pos),]
	
	negative$tag = rep(0,nrow(negative))
	positive$tag = rep(1,nrow(positive))
	
	data = rbind(positive,negative)
	
	if(i == 1)
	{
		row.names(res) = row.names(data)
		res[,i] = data$V4
	}else{
	
		res[,i] = data$V4
	}

}
raw.res = res
res = apply(res,2,function(x) as.numeric(as.vector(x)))
#res = sign(res)
tag = c(rep(1,2000),rep(0,2000))
data = cbind(tag,res)
data = as.data.frame(data)
data$tag = factor(data$tag, levels =c(0,1))
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_TF_integration.Rda"
save.image(file = myoutf)

#########################################################
#load Histone modification profile
#########################################################
rm(list=ls())
#load the tag
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_neg.Rda")

pos = validation_pos
neg = validation_neg

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_counts_bins_logFC/"
myinf1 =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/Histone_file_annotation_bigwig.txt"

files = list.files(mydir)
name = gsub(".Rda","",files)
name = gsub("_ENCFF172DEA_positive","",name)
name = gsub("_ENCFF172DEA_negative","",name)
name = unique(name)
info = read.table(myinf1,sep="\t",quote=NULL,stringsAsFactors=F)
info = info[intersect(name,row.names(info)),]
colnames(info)= c("target","file_accession")

library(randomForest)
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Multi_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Uni_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/moveme.R")

AUC_score = list()

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_positive.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_negative.Rda")



rm(positive)
rm(negative)

res = matrix(0, 4000, nrow(info))
colnames(res) = row.names(info)
res = as.data.frame(res)

for(i in 1 : nrow(info))
{
	cat("\r",i)
	
	tmpinf1 = paste0(mydir,row.names(info)[i],"_ENCFF172DEA_negative.Rda")
	load(tmpinf1)
	negative = histone
	negative = negative[row.names(neg),]
	
	tmpinf2 = paste0(mydir,row.names(info)[i],"_ENCFF172DEA_positive.Rda")
	load(tmpinf2)
	positive = histone
	positive = positive[row.names(pos),]
	
	negative$tag = rep(0,nrow(negative))
	positive$tag = rep(1,nrow(positive))
	
	data = rbind(positive,negative)
	
	if(i == 1)
	{
		row.names(res) = row.names(data)
		res[,i] = data$V4
	}else{
	
		res[,i] = data$V4
	}

}
raw.res = res
res = apply(res,2,function(x) as.numeric(as.vector(x)))
#res = sign(res)
tag = c(rep(1,2000),rep(0,2000))
data = cbind(tag,res)
data = as.data.frame(data)
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_Histone_integration.Rda"
save.image(file = myoutf)

#########################################################
#load HI-C profile
#########################################################
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Hi_C_contact_in_anchor/"
myinf1 =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/Hi_C_file_annotation_hic.txt"

files = list.files(mydir)
name = gsub("_contact_negative__normalized_counts.txt","",files)
name = gsub("_contact_positive__normalized_counts.txt","",name)
name = unique(name)
info = read.table(myinf1,sep="\t",quote=NULL,stringsAsFactors=F)
row.names(info) = info$file_accession
info = info[name,]
info = info[,c("target","dataset_description")]

library(randomForest)
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Multi_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Uni_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/moveme.R")

AUC_score = list()

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_positive.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_negative.Rda")

#load the tag
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_neg.Rda")

pos = validation_pos
neg = validation_neg

bin.size = 100

#begin to calcualte the pos
pos$V2 = as.numeric(as.vector(pos$V2))
pos$V3 = as.numeric(as.vector(pos$V3))
	
pos$sta = ceiling((pos$V2+1)/bin.size)
pos$end = ceiling((pos$V3+1)/bin.size)
	
pos$ID = paste0(pos$V1,"_",pos$sta)
row.names(pos) = pos$ID
	
#begin to calulate the neg
	
neg$V2 = as.numeric(as.vector(neg$V2))
neg$V3 = as.numeric(as.vector(neg$V3))
	
neg$sta = ceiling((neg$V2+1)/bin.size)
neg$end = ceiling((neg$V3+1)/bin.size)
	
neg$ID = paste0(neg$V1,"_",neg$sta)
row.names(neg) = neg$ID

rm(positive)
rm(negative)

res = matrix(0, 10000, nrow(info))
colnames(res) = row.names(info)
res = as.data.frame(res)

for(i in 1 : nrow(info))
{
	cat("\r",i)
	
	tmpinf1 = paste0(mydir,row.names(info)[i],"_contact_negative__normalized_counts.txt")
	negative = read.table(tmpinf1,sep="\t",quote=NULL)
	row.names(negative) = negative$ID
	negative = negative[row.names(neg),]
	
	tmpinf2 = paste0(mydir,row.names(info)[i],"_contact_positive__normalized_counts.txt")
	positive = read.table(tmpinf2,sep="\t",quote=NULL)
	row.names(positive) = positive$ID
	positive = positive[row.names(pos),]
	
	negative$tag = rep(0,nrow(negative))
	positive$tag = rep(1,nrow(positive))
	
	data = rbind(positive,negative)
	
	if(i == 1)
	{
		row.names(res) = row.names(data)
		res[,i] = data$value
	}else{
	
		res[,i] = data$value
	}

}
raw.res = res
res = apply(res,2,function(x) as.numeric(as.vector(x)))
#res = sign(res)
tag = c(rep(1,5000),rep(0,5000))
data = cbind(tag,res)
data = as.data.frame(data)
row.names(data) = c(row.names(pos),row.names(neg))
tag = is.na(data)
data[tag]=0
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_HIC_integration.Rda"
save.image(file = myoutf)

#########################################################
#load motif profile
#########################################################
rm(list=ls())
mydir = "/lorax/chenglab/cc59/PubDat/organisms/human/motif/hsa_done/"
files = list.files(mydir)

library(randomForest)
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Multi_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Uni_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/moveme.R")

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_neg.Rda")

bin.size = 100

pos = validation_pos
neg = validation_neg

#begin to calcualte the pos
pos$V2 = as.numeric(as.vector(pos$V2))
pos$V3 = as.numeric(as.vector(pos$V3))
	
pos$sta = ceiling((pos$V2+1)/bin.size)
pos$end = ceiling((pos$V3+1)/bin.size)
	
pos$ID = paste0(pos$V1,"_",pos$sta)
row.names(pos) = pos$ID
	
#begin to calulate the neg
	
neg$V2 = as.numeric(as.vector(neg$V2))
neg$V3 = as.numeric(as.vector(neg$V3))
	
neg$sta = ceiling((neg$V2+1)/bin.size)
neg$end = ceiling((neg$V3+1)/bin.size)
	
neg$ID = paste0(neg$V1,"_",neg$sta)
row.names(neg) = neg$ID
	
#create the matrix
motif_name = gsub(".txt","",files)

pos_count = matrix(0, nrow(pos), length(motif_name))	
row.names(pos_count) = row.names(pos)
colnames(pos_count) = motif_name
pos_count = as.data.frame(pos_count)

neg_count = matrix(0, nrow(neg), length(motif_name))
row.names(neg_count) = row.names(neg)
colnames(neg_count) = motif_name
neg_count = as.data.frame(neg_count)
	
for(i in 1 : length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir,files[i])
	tmp = read.table(tmpinf,sep="\t",quote=NULL)
	
	pos_com = intersect(row.names(pos_count),tmp$V1)
	neg_com = intersect(row.names(neg_count),tmp$V1)
	
	pos_count[pos_com,i]= 1
	neg_count[neg_com,i]= 1
	
}

tag = c(rep(1,nrow(pos_count)),rep(0,nrow(neg_count)))
res = rbind(pos_count, neg_count)
colnames(res) = paste0("F_",seq(1,ncol(pos_count)))

data = cbind(tag,res)
data = as.data.frame(data)
data$tag = factor(data$tag, levels =c(0,1))

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_motif_integration.Rda"
save.image(file = myoutf)
#########################################################
#load sequence profile
##########################################################
rm(list=ls())
library("Biostrings")
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Sequence_info/ENCFF172DEA_positive_1D_RF_score.fa"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Sequence_info/ENCFF172DEA_positive_1D_NN_score.fa"
myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Sequence_info/ENCFF172DEA_negative_1D_RF_score.fa"
myinf4 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Sequence_info/ENCFF172DEA_negative_1D_NN_score.fa"

###########################
#process the positive case
###########################

RF = readLines(myinf1)
tag2 = seq(2,10000,2)
tag1 = seq(1,10000,2)
RF_ID = RF[tag1]
RF_value = RF[tag2]
res1 = as.data.frame(cbind(RF_ID,RF_value))

RF = readLines(myinf2)
tag2 = seq(2,10000,2)
tag1 = seq(1,10000,2)
RF_ID = RF[tag1]
RF_value = RF[tag2]
res2 = as.data.frame(cbind(RF_ID,RF_value))

row.names(res1) = res1$RF_ID
row.names(res2) = res2$RF_ID

com = intersect(row.names(res1),row.names(res2))
res1 = res1[com,]
res2 = res2[com,]

pos = data.frame(ID = row.names(res1),RF_value = as.numeric(as.vector(res1$RF_value)), NN_value = as.numeric(as.vector(res2$RF_value)))
row.names(pos) = pos$ID
row.names(pos) = gsub(">","",row.names(pos))


##########################
#process the negative case
##########################
RF = readLines(myinf3)
tag2 = seq(2,10000,2)
tag1 = seq(1,10000,2)
RF_ID = RF[tag1]
RF_value = RF[tag2]
res1 = as.data.frame(cbind(RF_ID,RF_value))

RF = readLines(myinf4)
tag2 = seq(2,10000,2)
tag1 = seq(1,10000,2)
RF_ID = RF[tag1]
RF_value = RF[tag2]
res2 = as.data.frame(cbind(RF_ID,RF_value))

row.names(res1) = res1$RF_ID
row.names(res2) = res2$RF_ID

com = intersect(row.names(res1),row.names(res2))
res1 = res1[com,]
res2 = res2[com,]

neg = data.frame(ID = row.names(res1),RF_value = as.numeric(as.vector(res1$RF_value)), NN_value = as.numeric(as.vector(res2$RF_value)))
row.names(neg) = neg$ID
row.names(neg) = gsub(">","",row.names(neg))

res = rbind(pos,neg)
tag = c(rep(1,nrow(pos)),rep(0,nrow(neg)))
res = cbind(tag, res)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_sequence_information.txt"
write.table(res,myoutf,sep="\t",quote=F)

#########################################################
#load conservation file
##########################################################
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Sequence_info/validation_pos_conScore.txt"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Sequence_info/validation_neg_conScore.txt"

res1 = read.table(myinf1,quote=NULL,fill=T)
res2 = read.table(myinf2,quote=NULL,fill=T)

res = rbind(res1,res2)
row.names(res) = paste0(res[,1],"_",res[,2],"_",res[,3])

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_conScore_information.txt"
write.table(res,myoutf,sep="\t",quote=F)


############################
#Histone only
################################
rm(list=ls())

library(randomForest)

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_Histone_integration.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_TF_integration.Rda"
myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_motif_integration.Rda"

load(myinf1)
res1 = data
row.names(res1) = c(row.names(positive),row.names(negative))
res1 = res1[complete.cases(res1),]
res1$tag = factor(res1$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

fit = ROC(res1) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_integration_prediction.Rda"
save.image(file=myoutf)

############################
#TF only
################################
rm(list=ls())

library(randomForest)

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_Histone_integration.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_TF_integration.Rda"
myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_motif_integration.Rda"

load(myinf2)
res2 = data
row.names(res2) = c(row.names(positive),row.names(negative))
res2 = res2[complete.cases(res2),]
res2$tag = factor(res2$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

fit = ROC(res2) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda"
save.image(file=myoutf)

############################
#motif only
################################
rm(list=ls())

library(randomForest)

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_Histone_integration.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_TF_integration.Rda"
myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_motif_integration.Rda"

load(myinf3)
res3 = data
#row.names(res3) = c(row.names(positive),row.names(negative))
res3 = res3[complete.cases(res3),]
res3$tag = factor(res3$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

fit = ROC(res3) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/motif_integration_prediction.Rda"
save.image(file=myoutf)

############################
#HI-C only
################################
rm(list=ls())

library(randomForest)

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_Histone_integration.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_TF_integration.Rda"
myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_motif_integration.Rda"
myinf4 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_HIC_integration.Rda"

load(myinf4)
res4 = data
#row.names(res4) = c(row.names(positive),row.names(negative))
#res4 = res4[complete.cases(res4),]
res4$tag = factor(res4$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

fit = ROC(res4) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/HIC_integration_prediction.Rda"
save.image(file=myoutf)

############################
#Sequence only
################################
rm(list=ls())
library(randomForest)
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_sequence_information.txt"
data = read.table(myinf,sep="\t",quote=NULL)
data = data[,c("tag", "RF_value", "NN_value")]

data$tag = factor(data$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

fit = ROC(data) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Sequence_integration_prediction.Rda"
save.image(file=myoutf)

 1-fit[[3]]
 [1] 0.6050160 0.6041200 0.5975904 0.5986564 0.5990663 0.6021571 0.6010092
 [8] 0.6013708 0.6021426 0.6002387

############################
#Conservation only
################################
rm(list=ls())
library(randomForest)
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_conScore_information.txt"
data = read.table(myinf,sep="\t",quote=NULL)
tag = c(rep(1,5000),rep(0,5000))
data = cbind(tag,data)

data$tag = factor(data$tag, levels =c(0,1))
data = data[,c("tag", "conScore")]
data = data[complete.cases(data),]

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_Uni_variable_10_CV.R")

fit = ROC_uni(data) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/conScore_integration_prediction.Rda"
save.image(file=myoutf)
#######################################################################################################################
#Histone + TF
#######################################################################################################################
rm(list=ls())

library(randomForest)

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_Histone_integration.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_TF_integration.Rda"
myinf3 =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_conScore_information.txt"

load(myinf1)
res1 = data
row.names(res1) = c(row.names(positive),row.names(negative))
res1 = res1[complete.cases(res1),]

load(myinf2)
res2 = data
row.names(res2) = c(row.names(positive),row.names(negative))
res2 = res2[complete.cases(res2),]

com = Reduce(intersect,list(row.names(res1),row.names(res2)))
res1 = res1[com,]
res2 = res2[com,]

tag = res1[,1]
res1 = res1[,2:ncol(res1)]
res2 = res2[,2:ncol(res2)]

res = cbind(tag, res1,res2)
res$tag = factor(res$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

fit = ROC(res) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_TF_integration_prediction.Rda"
save.image(file=myoutf)

#######################################################################################################################
#Histone + TF + conScore
#######################################################################################################################
rm(list=ls())

library(randomForest)

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_Histone_integration.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_TF_integration.Rda"
myinf3 =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_conScore_information.txt"

load(myinf1)
res1 = data
row.names(res1) = c(row.names(positive),row.names(negative))
res1 = res1[complete.cases(res1),]

load(myinf2)
res2 = data
row.names(res2) = c(row.names(positive),row.names(negative))
res2 = res2[complete.cases(res2),]

#conservation score
res3 = read.table(myinf3,sep="\t",quote=NULL)
res3 = res3[complete.cases(res3),]

com = Reduce(intersect,list(row.names(res1),row.names(res2),row.names(res3)))
res1 = res1[com,]
res2 = res2[com,]
res3 = res3[com,]

tag = res1[,1]
res1 = res1[,2:ncol(res1)]
res2 = res2[,2:ncol(res2)]
conScore = res3[,"conScore"]

res = cbind(tag, res1,res2,conScore)
res$tag = factor(res$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

fit = ROC(res) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_TF_conScore_integration_prediction.Rda"
save.image(file=myoutf)

#######################################################################################################################
#Histone + TF + conScore + sequence
#######################################################################################################################
rm(list=ls())

library(randomForest)

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_Histone_integration.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_TF_integration.Rda"
myinf3 =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_conScore_information.txt"
myinf4 =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_sequence_information.txt"

load(myinf1)
res1 = data
row.names(res1) = c(row.names(positive),row.names(negative))
res1 = res1[complete.cases(res1),]

load(myinf2)
res2 = data
row.names(res2) = c(row.names(positive),row.names(negative))
res2 = res2[complete.cases(res2),]

#conservation score
res3 = read.table(myinf3,sep="\t",quote=NULL)
res3 = res3[complete.cases(res3),]

#sequence information
res4 = read.table(myinf4,sep="\t",quote=NULL)
res4 = res4[complete.cases(res4),]

com = Reduce(intersect,list(row.names(res1),row.names(res2),row.names(res3),row.names(res4)))
res1 = res1[com,]
res2 = res2[com,]
res3 = res3[com,]
res4 = res4[com,]

tag = res1[,1]
res1 = res1[,2:ncol(res1)]
res2 = res2[,2:ncol(res2)]
conScore = res3[,"conScore"]
res4 = res4[,c("RF_value", "NN_value")]

res = cbind(tag, res1,res2,conScore,res4)
res$tag = factor(res$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

fit = ROC(res) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_TF_conScore_sequence_integration_prediction.Rda"
save.image(file=myoutf)
#######################################################################################################################
rm(list=ls())

library(randomForest)

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_Histone_integration.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_TF_integration.Rda"
myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_motif_integration.Rda"
myinf4 =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_conScore_information.txt"

load(myinf1)
res1 = data
row.names(res1) = c(row.names(positive),row.names(negative))
res1 = res1[complete.cases(res1),]

load(myinf2)
res2 = data
row.names(res2) = c(row.names(positive),row.names(negative))
res2 = res2[complete.cases(res2),]

load(myinf3)
res3 = data
res3 = res3[complete.cases(res3),]

#process the row names of res3
tmp = row.names(res3)
chr = lapply(tmp, function(x) strsplit(x,"_")[[1]][1])
chr = as.vector(unlist(chr))

sta = lapply(tmp, function(x) strsplit(x,"_")[[1]][2])
sta = as.numeric(as.vector(unlist(sta)))
sta = (sta-1)*100

end = sta+99

my_bin_name = mapply(function(x,y,z) paste0(x,"_",y,"_",z), chr, sta, end)
row.names(res3) = my_bin_name

#conservation score
res4 = read.table(myinf4,sep="\t",quote=NULL)
res4 = res4[complete.cases(res4),]

com = Reduce(intersect,list(row.names(res1),row.names(res2),row.names(res3),row.names(res4)))
res1 = res1[com,]
res2 = res2[com,]
res3 = res3[com,]
res4 = res4[com,]

tag = res1[,1]
res1 = res1[,2:ncol(res1)]
res2 = res2[,2:ncol(res2)]
res3 = res3[,2:ncol(res3)]
conScore = res4[,"conScore"]

res = cbind(tag, res1,res2,res3,conScore)
res$tag = factor(res$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

fit = ROC(res) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Without_seq_integration_prediction.Rda"
save.image(file=myoutf)


#######################################################################################################################
#integration with sequence
#######################################################################################################################
rm(list=ls())

library(randomForest)

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_Histone_integration.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_TF_integration.Rda"
myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_motif_integration.Rda"
myinf4 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_sequence_information.txt"
myinf5 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_conScore_information.txt"

load(myinf1)
res1 = data
row.names(res1) = c(row.names(positive),row.names(negative))
res1 = res1[complete.cases(res1),]

load(myinf2)
res2 = data
row.names(res2) = c(row.names(positive),row.names(negative))
res2 = res2[complete.cases(res2),]

load(myinf3)
res3 = data
res3 = res3[complete.cases(res3),]

#process the row names of res3
tmp = row.names(res3)
chr = lapply(tmp, function(x) strsplit(x,"_")[[1]][1])
chr = as.vector(unlist(chr))

sta = lapply(tmp, function(x) strsplit(x,"_")[[1]][2])
sta = as.numeric(as.vector(unlist(sta)))
sta = (sta-1)*100

end = sta+99

my_bin_name = mapply(function(x,y,z) paste0(x,"_",y,"_",z), chr, sta, end)

row.names(res3) = my_bin_name

res4 = read.table(myinf4,sep="\t",quote=NULL)
res4 = res4[complete.cases(res4),]

res5 = read.table(myinf5,sep="\t",quote=NULL)
res5 = res5[complete.cases(res5),]

com = Reduce(intersect,list(row.names(res1),row.names(res2),row.names(res3),row.names(res4), row.names(res5)))
res1 = res1[com,]
res2 = res2[com,]
res3 = res3[com,]
res4 = res4[com,]
res5 = res5[com,]

tag = res1[,1]
res1 = res1[,2:ncol(res1)]
res2 = res2[,2:ncol(res2)]
res3 = res3[,2:ncol(res3)]
res4 = res4[,c("RF_value","NN_value")]
conScore = res5[,"conScore"]

res = cbind(tag, res1,res2,res3,res4,conScore)
res$tag = factor(res$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

fit = ROC(res) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Final_integration_prediction.Rda"
save.image(file=myoutf)

##################################################
#Promoter vs non-promoter
#################################################
rm(list=ls())
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Final_integration_prediction.Rda")
raw.data = res

#load the tag
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda")
pos = validation_pos
#bin.size = 100
#pos$V2 = as.numeric(as.vector(pos$V2))
#pos$sta = ceiling((pos$V2+1)/bin.size)
#row.names(pos) = paste0(pos$V1,"_",pos$sta)

data = raw.data
com = intersect(row.names(pos),row.names(data))
data = data[com,]

tmpinf =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_positive_annotation.txt"
info = read.table(tmpinf,sep="\t",quote=NULL)
#info$sta = ceiling((info$start+1)/bin.size)
row.names(info) = paste0(info$seqnames,"_",info$start,"_",info$end)

com = intersect(row.names(info),row.names(data))
data= data[com,]
info = info[com,]

#get the positive data
tag = data$tag == "1"
data = data[tag,]
data$tag = as.vector(data$tag)
data$tag = rep(0,nrow(data))

tag = grep("Promoter",info$annotation)
data[tag,"tag"] = 1
data$tag = factor(data$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

fit = ROC(data) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Final_integration_prediction_promoter_distal.Rda"
save.image(file=myoutf)

######################################################
#Leave each unit out for testing Histone
#####################################################
#submit to server

rm(list=ls())
library(randomForest)

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_integration_prediction.Rda")

#relative importance
res = fit[[6]]
res = apply(res,1,mean)
res = res[order(res,decreasing=T)]
tar = names(res)

for(i in 1 : (length(tar)-1))
{
	cat("\r",i)
	
	sam_for_taken = tar[1:i]
	tag = which(colnames(data) %in% sam_for_taken)
	test = data[,-tag]
	
	test  = test[complete.cases(test),]
	test[,"tag"] = factor(test[,"tag"], levels =c(0,1))

	if(i==38)
	{
		source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_Uni_variable_10_CV.R")
		fit = ROC_uni(test)
		myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_KO_RF/GM12878_Order_",i,"_",tar[i],"_KO_histone.Rda")
		save(fit, file = myoutf)
		
	}else{
	
		source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")	
		fit = ROC(test) 
		rank = fit[[6]]
		rank = apply(rank,1,mean)
		rank = rank[order(rank,decreasing=T)]
		tar[i+1] = names(rank)[1]	
		
		myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_KO_RF/GM12878_Order_",i,"_",tar[i],"_KO_histone.Rda")
		save(fit, file = myoutf)
		}
	
	
	#begin to calculate the rank
}

######################################################
#Leave each unit out for testing TF
#####################################################
#submit to server

#PBS -l walltime=200:00:00 -l feature=celln
/usr/bin/R --vanilla <<EOF

rm(list=ls())
library(randomForest)

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda")

#relative importance
res = fit[[6]]
res = apply(res,1,mean)
res = res[order(res,decreasing=T)]
tar = names(res)

for(i in 1 : (length(tar)-1))
{
	cat("\r",i)
	
	sam_for_taken = tar[1:i]
	tag = which(colnames(data) %in% sam_for_taken)
	test = data[,-tag]
	
	test  = test[complete.cases(test),]
	test[,"tag"] = factor(test[,"tag"], levels =c(0,1))
	
	if(i==466)
	{
		source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_Uni_variable_10_CV.R")
		fit = ROC_uni(test)
		myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_KO_RF/GM12878_Order_",i,"_",tar[i],"_KO_TF.Rda")
		save(fit, file = myoutf)
	}else{
	
		source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")	
		fit = ROC(test) 
		
		myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_KO_RF/GM12878_Order_",i,"_",tar[i],"_KO_TF.Rda")
		save(fit, file = myoutf)
	
		#begin to calculate the rank
	
		rank = fit[[6]]
		rank = apply(rank,1,mean)
		rank = rank[order(rank,decreasing=T)]
		tar[i+1] = names(rank)[1]	
	}
	
	
}

EOF







