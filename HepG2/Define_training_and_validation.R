[1]Merge all the bins of histones
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Histone_counts_bins_logFC/"
files = list.files(mydir)
tag = grep("positive",files)
files = files[tag]

for(i in 1: length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir, files[i])
	load(tmpinf)
	
	if(i == 1)
	{
		bins = histone$ID
	}else{
	
		bins = intersect(bins, histone$ID)
	}

}

length(bins)
[1] 882920

positive_bins = bins
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Total_bins/Pos_bin_Histone.Rda"
save(positive_bins,file = myoutf)


#begin to test the shared number (Positive)
#begin to test the shared number (Positive)
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Histone_counts_bins_logFC/"
files = list.files(mydir)
tag = grep("negative",files)
files = files[tag]

for(i in 1: length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir, files[i])
	load(tmpinf)
	
	if(i == 1)
	{
		bins = histone$ID
	}else{
	
		bins = intersect(bins, histone$ID)
	}

}
length(bins)
[1] 882884

negative_bins = bins
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Total_bins/Neg_bin_Histone.Rda"
save(negative_bins,file = myoutf)

[2]Merge all the bins of TFs
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_counts_bins_logFC/"
files = list.files(mydir)
tag = grep("positive",files)
files = files[tag]

for(i in 1: length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir, files[i])
	load(tmpinf)
	
	if(i == 1)
	{
		bins = histone$ID
	}else{
	
		bins = intersect(bins, histone$ID)
	}
	
	if(length(bins) < 1)
	{
		command = paste0("Warning job",i)
		print(command)
	}
}
length(bins)
[1]  873810

positive_bins = bins
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Total_bins/Pos_bin_TF.Rda"
save(positive_bins,file = myoutf)

#begin to test the shared number (Positive)
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_counts_bins_logFC/"
files = list.files(mydir)
tag = grep("negative",files)
files = files[tag]

for(i in 1: length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir, files[i])
	load(tmpinf)
	
	if(i == 1)
	{
		bins = histone$ID
	}else{
	
		bins = intersect(bins, histone$ID)
	}

}

length(bins)
[1] 895704
negative_bins = bins
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Total_bins/Neg_bin_TF.Rda"
save(negative_bins,file = myoutf)

[3]Beign to calculate the all bins in Histone
rm(list=ls())
#load the tag
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Histone_counts_bins_logFC/"
myinf1 =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/Histone_file_annotation_bigwig.txt"

files = list.files(mydir)
name = gsub(".Rda","",files)
name = gsub("_ENCFF356TXH_positive","",name)
name = gsub("_ENCFF356TXH_negative","",name)
name = unique(name)
info = read.table(myinf1,sep="\t",quote=NULL,stringsAsFactors=F)
info = info[intersect(name,row.names(info)),]
colnames(info)= c("target","file_accession")

library(randomForest)
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Multi_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Uni_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/moveme.R")

AUC_score = list()

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Total_bins/Pos_bin_Histone.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Total_bins/Neg_bin_Histone.Rda")

res = matrix(0, 1765804, nrow(info))
colnames(res) = row.names(info)
res = as.data.frame(res)

for(i in 1 : nrow(info))
{
	cat("\r",i)
	
	tmpinf1 = paste0(mydir,row.names(info)[i],"_ENCFF356TXH_negative.Rda")
	load(tmpinf1)
	negative = histone
	negative = negative[negative_bins,]
	
	tmpinf2 = paste0(mydir,row.names(info)[i],"_ENCFF356TXH_positive.Rda")
	load(tmpinf2)
	positive = histone
	positive = positive[positive_bins,]
	
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
tag = c(rep(1,length(positive_bins)),rep(0,length(negative_bins)))
data = cbind(tag,res)
data = as.data.frame(data)
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_Histone_integration_Total.Rda"
save.image(file = myoutf)

[4]Beign to calculate the all bins in TF
rm(list=ls())
#load the tag
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_counts_bins_logFC/"
myinf1 =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/TF_file_annotation_bigWig.txt"
 
files = list.files(mydir)
name = gsub(".Rda","",files)
name = gsub("_ENCFF356TXH_positive","",name)
name = gsub("_ENCFF356TXH_negative","",name)
name = unique(name)
info = read.table(myinf1,sep="\t",quote=NULL,stringsAsFactors=F)
info = info[intersect(name,row.names(info)),]
colnames(info)= c("target","file_accession")

library(randomForest)
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Multi_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Uni_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/moveme.R")

AUC_score = list()

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Total_bins/Pos_bin_TF.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Total_bins/Neg_bin_TF.Rda")

res = matrix(0, 1766007, nrow(info))
colnames(res) = row.names(info)
res = as.data.frame(res)

for(i in 1 : nrow(info))
{
	cat("\r",i)
	
	tmpinf1 = paste0(mydir,row.names(info)[i],"_ENCFF356TXH_negative.Rda")
	load(tmpinf1)
	negative = histone
	negative = negative[negative_bins,]
	
	tmpinf2 = paste0(mydir,row.names(info)[i],"_ENCFF356TXH_positive.Rda")
	load(tmpinf2)
	positive = histone
	positive = positive[positive_bins,]
	
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
tag = c(rep(1,length(positive_bins)),rep(0,length(negative_bins)))
#data = cbind(tag,res)
#data = as.data.frame(data)

res = cbind(tag,res)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_TF_integration_Total.Rda"
save.image(file = myoutf)


myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_TF_integration_Total_data_frame_only.Rda"
save(res,file = myoutf)

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_TF_integration_Total_data_frame_only.Rda")
data = res
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_TF_integration_Total_data_frame_only.Rda"
save(data,file = myoutf)

[5]Define the shared region between TF and Histone
rm(list=ls())
myinf1= "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_Histone_integration_Total.Rda"
#myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_TF_integration_Total.Rda"
myinf2 =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_TF_integration_Total_data_frame_only.Rda"

load(myinf1)
row.names(data) = row.names(raw.res)
res1 = data

load(myinf2)
#row.names(data) = row.names(raw.res)
res2 = data

com = intersect(row.names(res1),row.names(res2))
res1 = res1[com,]
res2 = res2[com,]

myoutf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_Histone_integration_Core.Rda"
save(res1, file = myoutf1)

myoutf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_TF_integration_Core.Rda"
save(res2, file = myoutf2)


[1]sequence define positive sequence #Use current histone as the golden standard for testing since its high AUC
rm(list=ls())
library("Biostrings")
library("reticulate")
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_bin_sequence/ENCFF356TXH_positive_peak.fa"
fastaFile <- readDNAStringSet(myinf)

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_Histone_integration_Core.Rda")
tag = res1$tag == "1"
positive = res1[tag,]

#take the shared bins and fasta files nmaes
com =intersect(row.names(positive), names(fastaFile))
positive = positive[com,]
fastaFile = fastaFile[com]
raw.positive = positive

tag = sample(length(com),2000)
validation_pos_tag = row.names(positive)[tag]

#restrain the shared positive
positive = raw.positive
ID = row.names(positive)
com = intersect(ID, names(fastaFile))
positive = positive[com,]
fastaFile = fastaFile[com]
ID = com

#pick up the label
total_seq = ID
training_pos_tag = setdiff(total_seq, validation_pos_tag)

#begin to get the sequence file
tag = sample(seq(1,length(training_pos_tag)),20000)

ID_for_training = training_pos_tag[tag]

training_fasta = fastaFile[ID_for_training]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_bin_sequence/ENCFF356TXH_positive_training.fa"
writeXStringSet(training_fasta, myoutf)

ID_for_validation = validation_pos_tag
validation_fasta = fastaFile[ID_for_validation]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_bin_sequence/ENCFF356TXH_positive_validation.fa"
writeXStringSet(validation_fasta, myoutf)

#begin to get the bed files

training_pos = positive[ID_for_training,]
validation_pos = positive[ID_for_validation,]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/training_pos_bed_file.txt"
write.table(training_pos, myoutf, col.names=F,row.names=F,quote=F)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/validation_pos_bed_file.txt"
write.table(validation_pos, myoutf, col.names=F,row.names=F,quote=F)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/HepG2_training_and_validation_pos.Rda"
save.image(file = myoutf)

[3]sequence define negative sequence
rm(list=ls())
library("Biostrings")
library("reticulate")
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_bin_sequence/ENCFF356TXH_negative_peak.fa"
fastaFile <- readDNAStringSet(myinf)

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_Histone_integration_Core.Rda")
tag = res1$tag == "0"
negative = res1[tag,]

#take the shared bins and fasta files nmaes
com =intersect(row.names(negative), names(fastaFile))
negative = negative[com,]
fastaFile = fastaFile[com]
raw.negative = negative

tag = sample(length(com),2000)
validation_neg_tag = row.names(negative)[tag]

#restrain the peak that can be plot the sequence
negative = raw.negative
ID = row.names(negative)
com = intersect(ID, names(fastaFile))
negative = negative[com,]
fastaFile = fastaFile[com]
ID = com

#pick up the label
total_seq = ID
training_neg_tag = setdiff(total_seq, validation_neg_tag)

#begin to get sequence file
tag = sample(seq(1,length(training_neg_tag)),20000)

ID_for_training = training_neg_tag[tag]

training_fasta = fastaFile[ID_for_training]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_bin_sequence/ENCFF356TXH_negative_training.fa"
writeXStringSet(training_fasta, myoutf)

ID_for_validation = validation_neg_tag
validation_fasta = fastaFile[ID_for_validation]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_bin_sequence/ENCFF356TXH_negative_validation.fa"
writeXStringSet(validation_fasta, myoutf)

#begin to get the bed files
training_neg = negative[ID_for_training,]
validation_neg = negative[ID_for_validation,]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/training_neg_bed_file.txt"
write.table(training_neg, myoutf, col.names=F,row.names=F,quote=F)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/validation_neg_bed_file.txt"
write.table(validation_neg, myoutf, col.names=F,row.names=F,quote=F)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/HepG2_training_and_validation_neg.Rda"
save.image(file = myoutf)


[4]bed profile prepration

cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/

cut -f 1,2,3 training_pos_bed_file.txt > training_pos.bed
cut -f 1,2,3 validation_pos_bed_file.txt > validation_pos.bed

cut -f 1,2,3 training_neg_bed_file.txt > training_neg.bed
cut -f 1,2,3 validation_neg_bed_file.txt > validation_neg.bed

