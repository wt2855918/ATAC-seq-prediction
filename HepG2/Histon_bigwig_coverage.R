[1]Begin to calculate the coverage 
rm(list=ls())
work_mydir  = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_Histone_from_bigWig/"
dir.create(work_mydir)
setwd(work_mydir)

sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/Histone_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files,fix=T)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

myoutdir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Histone_logFC/"
dir.create(myoutdir)

for(k in 1 : length(bam_files))
{
	cat("\r",k)
	myoutf1 = paste("job", k, ".sp", sep="")
	conOut = file(myoutf1, "w")
	curLine = c("#PBS -l walltime=10:00:00 -l vmem=16gb", 
				"", 
				"module load python/2.7-Anaconda", 
				"source activate SICER")
	writeLines(curLine, conOut)
	
	curLine = c("cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/Histone_ChIP_seq_bigWig/")
	writeLines(curLine, conOut)
	
	input_file = bam_files[k]
	
	output1 = paste0(myoutdir,sam[k],".npz")
	output2 = paste0(myoutdir,sam[k],".tab")
	
	curLine = paste0("multiBigwigSummary bins -bs 100 -b ", input_file," -out ", output1, " --outRawCounts ",output2)
	
	writeLines(curLine, conOut)	
	close(conOut)
}

#begin to submit job
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_Histone_from_bigWig/"
files = list.files(myDir1)

for(i in 1 : length(files))
{
	cat("\r",i)
	myinf = paste0(myDir1,files[i])
	command = paste0("qsub ", myinf)
	system(command)
}

##begin to check
#library(rtracklayer)
#myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/Histone_ChIP_seq_bigWig/ENCFF091LGA.bigWig"
#res = import(myinf, format = "bigWig", as = "GRanges")

[2]Begin to calculate shared bins template code for positive
rm(list=ls())
#define the function
f1 <- function(x, chr, sta, end) {paste(x[chr], x[sta]:x[end], sep="_")} #separate each base pair
f2 <- function(x, sig, sta, end) 
	{
		time = x[end]-x[sta]+1
		rep(as.character(x[sig]), time)
	}
	
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Histone_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]


myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_peak_processed_bed/ENCFF356TXH_positive.Rda"
load(myinf)
bin.size = 100

mytardir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Histone_counts_bins_logFC/"
my_nam = gsub(".tab","",files)

i= mykkk
cat("\r",i)
tmpinf = paste0(mydir,files[i])
data = read.table(tmpinf,sep="\t",quote=NULL)
data[,"V3"] = data[,"V3"] - 1

#delete the uncoverage region (due to hg19 version difference)
tag = is.na(data[,"V4"])
data = data[tag==0,]
	
tar = paste0("chr",c(seq(1,22,1),"X","Y"))
tag = which(data[,"V1"] %in% tar)
data = data[tag,]
	
histone = data
histone[,"ID"] = paste0(histone[,"V1"],"_",histone[,"V2"],"_",histone[,"V3"])

positive[,"ID"] = paste0(positive[,"V1"],"_",positive[,"V2"],"_",positive[,"V3"])

row.names(histone) = histone[,"ID"]
row.names(positive) = positive[,"ID"]
	
com = intersect(row.names(histone),row.names(positive))
	
histone = histone[com,]
	
myoutf = paste0(mytardir,my_nam[i],"_ENCFF356TXH_positive.Rda")
save(histone,file = myoutf)

[3]Begin to calculate the signal pos
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Histone_pos_ATAC/"
dir.create(myDir1)
mytemplate = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Template/HepG2/Template_for_counting_HepG2_Histone_sig_pos.R"
mysub = paste(myDir1, "submit.sp", sep="")
setwd(myDir1)

conIn = file(mytemplate, "r")
rawdata = readLines(conIn)
close(conIn)

fnum = 52
for(k in 1:fnum)
{
	data = rawdata
	se = grep("i= mykkk", data)
	tmp = data[se]
	tmp = gsub("mykkk", k, tmp)
	data[se] = tmp
	myoutf1 = paste("job", k, ".sp", sep="")
	conOut = file(myoutf1, "w")
	writeLines(data, conOut)
	close(conOut)
}

[4]Begin to calculate shared bins template code for negative
rm(list=ls())
#define the function
f1 <- function(x, chr, sta, end) {paste(x[chr], x[sta]:x[end], sep="_")} #separate each base pair
f2 <- function(x, sig, sta, end) 
	{
		time = x[end]-x[sta]+1
		rep(as.character(x[sig]), time)
	}
	
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Histone_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]


myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_peak_processed_bed/ENCFF356TXH_negative.Rda"
load(myinf)
bin.size = 100

mytardir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Histone_counts_bins_logFC/"
my_nam = gsub(".tab","",files)

i= mykkk
cat("\r",i)
tmpinf = paste0(mydir,files[i])
data = read.table(tmpinf,sep="\t",quote=NULL)
data[,"V3"] = data[,"V3"] - 1

#delete the uncoverage region (due to hg19 version difference)
tag = is.na(data[,"V4"])
data = data[tag==0,]
	
tar = paste0("chr",c(seq(1,22,1),"X","Y"))
tag = which(data[,"V1"] %in% tar)
data = data[tag,]
	
histone = data
histone[,"ID"] = paste0(histone[,"V1"],"_",histone[,"V2"],"_",histone[,"V3"])

negative[,"ID"] = paste0(negative[,"V1"],"_",negative[,"V2"],"_",negative[,"V3"])

row.names(histone) = histone[,"ID"]
row.names(negative) = negative[,"ID"]
	
com = intersect(row.names(histone),row.names(negative))
	
histone = histone[com,]
	
myoutf = paste0(mytardir,my_nam[i],"_ENCFF356TXH_negative.Rda")
save(histone,file = myoutf)

[5]Begin to calcualte the signal for negative
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Histone_neg_ATAC/"
dir.create(myDir1)
mytemplate = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Template/HepG2/Template_for_counting_HepG2_Histone_sig_neg.R"
mysub = paste(myDir1, "submit.sp", sep="")
setwd(myDir1)

conIn = file(mytemplate, "r")
rawdata = readLines(conIn)
close(conIn)

fnum = 52
for(k in 1:fnum)
{
	data = rawdata
	se = grep("i= mykkk", data)
	tmp = data[se]
	tmp = gsub("mykkk", k, tmp)
	data[se] = tmp
	myoutf1 = paste("job", k, ".sp", sep="")
	conOut = file(myoutf1, "w")
	writeLines(data, conOut)
	close(conOut)
}

[6]Model fiting
#begin to test the shared number (Positive)
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
#begin to test the shared number (Negative)
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

rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Histone_counts_bins_logFC/"
myinf1 =   "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/Histone_file_annotation_bigWig.txt"

files = list.files(mydir)
name = gsub(".Rda","",files)
name = gsub("_ENCFF356TXH_positive","",name)
name = gsub("_ENCFF356TXH_negative","",name)
name = unique(name)
info = read.table(myinf1,sep="\t",quote=NULL,stringsAsFactors=F)

library(randomForest)
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Multi_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Uni_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/moveme.R")

AUC_score = list()

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_peak_processed_bed/ENCFF356TXH_positive.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_peak_processed_bed/ENCFF356TXH_negative.Rda")


tag_pos = sample(seq(1,882920),2000)
tag_neg = sample(seq(1,882884),2000)

rm(positive)
rm(negative)

res = matrix(0, 4000, nrow(info))
colnames(res) = row.names(info)
res = as.data.frame(res)

for(i in 1 : nrow(info))
{
	cat("\r",i)
	
	tmpinf1 = paste0(mydir,row.names(info)[i],"_ENCFF356TXH_negative.Rda")
	load(tmpinf1)
	negative = histone
	negative = negative[tag_neg,]
	
	tmpinf2 = paste0(mydir,row.names(info)[i],"_ENCFF356TXH_positive.Rda")
	load(tmpinf2)
	positive = histone
	positive = positive[tag_pos,]
	
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

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(data)


1-fit[[3]]
 [1] 0.8815875 0.8722594 0.8824681 0.8742797 0.8712110 0.8721056 0.8675796
 [8] 0.8650334 0.8672085 0.8699019
  
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/HepG2_Histone_integration_AUC.Rda"
save.image(file = myoutf)

rank = fit[[6]]
rank = apply(rank,1,mean)
rank = rank[order(rank,decreasing=T)]
rank_info = info[names(rank),]




