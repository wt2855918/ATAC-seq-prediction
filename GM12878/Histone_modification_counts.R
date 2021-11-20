[1]Separate the bed graph files into bins 
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_positive.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_negative.Rda"

outdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/Bed_files_ATAC_seq/Split_bed/")
dir.create(outdir)

load(myinf1)

positive$V1 = as.character(as.vector(positive$V1))
positive$V2 = as.numeric(as.vector(positive$V2))
positive$V3 = as.numeric(as.vector(positive$V3))
positive = positive[order(positive$V2),]

mybin = ceiling(nrow(positive)/1000)
for(i in 1 : mybin)
{
	cat("\r",i)
	sta = 1000*(i-1) + 1
	end = 1000*i
	
	if(i ==  mybin )
	{
		end = nrow(positive)
	}
	
	res = positive[sta:end,]
	
	myoutf = paste0(outdir,"Positive_Bin_",i,".txt")
	write.table(res, myoutf,sep="\t",quote=F, row.names=F,col.names=F)

}
#####################################
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_negative.Rda"

outdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/Bed_files_ATAC_seq/Split_bed/")
dir.create(outdir)

load(myinf1)

negative$V1 = as.character(as.vector(negative$V1))
negative$V2 = as.numeric(as.vector(negative$V2))
negative$V3 = as.numeric(as.vector(negative$V3))
negative = negative[order(negative$V2),]

mybin = ceiling(nrow(negative)/1000)
for(i in 1 : mybin)
{
	cat("\r",i)
	sta = 1000*(i-1) + 1
	end = 1000*i
	
	if(i ==  mybin )
	{
		end = nrow(negative)
	}
	
	res = negative[sta:end,]
	
	myoutf = paste0(outdir,"Negative_Bin_",i,".txt")
	write.table(res, myoutf,sep="\t",quote=F, row.names=F,col.names=F)

}
#####################################
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Bed_files_ATAC_seq/Split_bed/"
files = list.files(mydir)
setwd("/lorax/chenglab/yanding/ATAC_seq_integration/Data/Bed_files_ATAC_seq/Split_bed/")
my_nam = gsub(".txt","",files)

for(i in 1 : length(files))
{
	cat("\r",i)
	input = paste0(my_nam[i],".txt")
	output = paste0(my_nam[i],".bed")
	command = paste0("cut -f 1,2,3 ",input," > ", output)
	system(command)
}



positive = positive[4001:5000,]
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Bed_files_ATAC_seq/ENCFF172DEA_positive.txt"
write.table(positive, myoutf,sep="\t",quote=F, row.names=F,col.names=F)

load(myinf2)
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Bed_files_ATAC_seq/ENCFF172DEA_negative.txt"
write.table(negative, myoutf,sep="\t",quote=F, row.names=F,col.names=F)

cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Bed_files_ATAC_seq/
cut -f 1,2,3 ENCFF172DEA_negative.txt > ENCFF172DEA_negative.bed
cut -f 1,2,3 ENCFF172DEA_positive.txt > ENCFF172DEA_positive.bed

[3]begin to quantify the counts
cd /lorax/chenglab/yanding/CoLab/SWI_SNF/A1a_phf10_180_AP1_sorted_bam/
multiBamSummary BED-file --BED /lorax/chenglab/yanding/CoLab/SWI_SNF/HiChIpper_res/Hi_Chip.bed --bamfiles D0_180_sorted.bam D0_A1a_sorted.bam D0_BRD9_sorted.bam D0_JunD_sorted.bam D0_pcJun_sorted.bam D20_155_sorted.bam D20_JunB_sorted.bam --labels D0_180 D0_A1a D0_BRD9 D0_JunD D0_pcJun D20_155 D20_JunB -o /lorax/chenglab/yanding/CoLab/SWI_SNF/HiChIP_Data/BAM_log2/Day_0_SWI_SNF_counts_summary.npz --outRawCounts /lorax/chenglab/yanding/CoLab/SWI_SNF/HiChIP_Data/BAM_log2/Day_0_SWI_SNF_counts_readCounts.tab
multiBamSummary BED-file --BED /lorax/chenglab/yanding/ATAC_seq_integration/Data/Bed_files_ATAC_seq/ENCFF172DEA_positive.bed --bamfiles D0_180_sorted.bam -o /lorax/chenglab/yanding/CoLab/SWI_SNF/HiChIP_Data/BAM_log2/test_summary.npz --outRawCounts /lorax/chenglab/yanding/CoLab/SWI_SNF/HiChIP_Data/BAM_log2/test.tab


cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/Histone_ChIP_seq/
multiBamSummary BED-file --BED /lorax/chenglab/yanding/ATAC_seq_integration/Data/Bed_files_ATAC_seq/ENCFF172DEA_positive.bed --bamfiles ENCFF958QVX.bam ENCFF019VEK.bam -o /lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_counts/ENCFF958QVX_ENCFF172DEA_bed_positive_summary.npz --outRawCounts /lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_counts/ENCFF958QVX_ENCFF172DEA_bed_positive.tab
bamCoverage --bam ENCFF958QVX.bam -o test.bedgraph



multiBamSummary BED-file --BED /lorax/chenglab/yanding/CoLab/SWI_SNF/HiChIpper_res/Hi_Chip.bed --bamfiles ENCFF958QVX.bam ENCFF019VEK.bam -o /lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_counts/ENCFF958QVX_ENCFF172DEA_bed_positive_summary.npz --outRawCounts /lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_counts/ENCFF958QVX_ENCFF172DEA_bed_positive.tab


cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/Histone_ChIP_seq/
"for i in ./*
do
samtools index $i
done"
module load python/2.7-Anaconda
source activate SICER

[2]Begin to calculate
rm(list=ls())
work_mydir  = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_Histone/"
dir.create(work_mydir)
setwd(work_mydir)

mydir  = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Bed_files_ATAC_seq/Split_bed/"
files = list.files(mydir)
tag = grep(".bed",files)
files = files[tag]
sam = gsub(".bed","",files)

sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/Histone_ChIP_seq/"
bam_files = list.files(sourcedir)
tag = grep("bai",bam_files)
bam_files = bam_files[-tag]
bam_files = paste(bam_files,collapse =" ")

myoutdir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_counts/"

for(k in 1 : length(files))
{
	cat("\r",k)
	myoutf1 = paste("job", k, ".sp", sep="")
	conOut = file(myoutf1, "w")
	curLine = c("#PBS -l walltime=10:00:00", 
				"", 
				"module load python/2.7-Anaconda", 
				"source activate SICER")
	writeLines(curLine, conOut)
	
	curLine = c("cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/Histone_ChIP_seq/")
	writeLines(curLine, conOut)
	
	bed_file = paste0(mydir,files[k])
	
	output1 = paste0(myoutdir,sam[k],"_summary.npz")
	output2 = paste0(myoutdir,sam[k],"_readcounts.tab")
	
	curLine = paste0("multiBamSummary BED-file --BED ", bed_file," --bamfiles ", bam_files, " -o ",output1," --outRawCounts ",output2)
	
	#curLine = paste("/opt/bowtie/2.2.7/bowtie2 --very-sensitive -x /lorax/chenglab/cc59/SofWar/reference/bowtie/hg19/hg19 -1 ", file1, " -2 ",file2," -S ",output)
	

	writeLines(curLine, conOut)	
	close(conOut)
	
	
}

#begin to submit job
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_Histone/"
files = list.files(myDir1)

for(i in 1 : length(files))
{
	cat("\r",i)
	myinf = paste0(myDir1,files[i])
	command = paste0("qsub ", myinf)
	system(command)
}

[3]Begin to calculate the full coverage
rm(list=ls())
work_mydir  = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_Histone/"
dir.create(work_mydir)
setwd(work_mydir)

sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/Histone_ChIP_seq/"
bam_files = list.files(sourcedir)
tag = grep("bai",bam_files)
bam_files = bam_files[-tag]
sam = gsub(".bam","",bam_files)

myoutdir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_counts/"

for(k in 1 : length(bam_files))
{
	cat("\r",k)
	myoutf1 = paste("job", k, ".sp", sep="")
	conOut = file(myoutf1, "w")
	curLine = c("#PBS -l walltime=10:00:00", 
				"", 
				"module load python/2.7-Anaconda", 
				"source activate SICER")
	writeLines(curLine, conOut)
	
	curLine = c("cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/Histone_ChIP_seq/")
	writeLines(curLine, conOut)
	
	input_file = bam_files[k]
	
	output1 = paste0(myoutdir,sam[k],".bedgraph")
	#output2 = paste0(myoutdir,sam[k],".bedgraph")
	
	curLine = paste0("bamCoverage --bam ", input_file," -o ", output1, " -of bedgraph --binSize 100 --effectiveGenomeSize 2864785220")
	
	#curLine = paste("/opt/bowtie/2.2.7/bowtie2 --very-sensitive -x /lorax/chenglab/cc59/SofWar/reference/bowtie/hg19/hg19 -1 ", file1, " -2 ",file2," -S ",output)
	

	writeLines(curLine, conOut)	
	close(conOut)
	
	
}

#begin to submit job
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_Histone/"
files = list.files(myDir1)

for(i in 1 : length(files))
{
	cat("\r",i)
	myinf = paste0(myDir1,files[i])
	command = paste0("qsub ", myinf)
	system(command)
}

[3]Begin to test the positive counts
rm(list=ls())
#define the function
f1 <- function(x, chr, sta, end) {paste(x[chr], x[sta]:x[end], sep="_")} #separate each base pair
f2 <- function(x, sig, sta, end) 
	{
		time = x[end]-x[sta]+1
		rep(as.character(x[sig]), time)
	}
	
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_counts/"
files = list.files(mydir)

myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_positive.Rda"
load(myinf)
bin.size = 100

mytardir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_counts_bins/"
my_nam = gsub(".bedgraph","",files)

for(i in 10 : length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir,files[i])
	data = read.table(tmpinf,sep="\t",quote=NULL)
	data$V2 = data$V2 +1 
	data$V3 = data$V3 +1 
	data$width =data$V3 - data$V2
	tag = data$width >= 100
	data = data[tag,]
	
	tar = paste0("chr",c(seq(1,22,1),"X","Y"))
	tag = which(data$V1 %in% tar)
	data = data[tag,]
	
	data$sta = ceiling((data$V2)/bin.size)
	data$end = floor((data$V3)/bin.size)
	
	tmp = unlist(apply(data, 1, f1, chr="V1", sta="sta", end="end"))
	tmp = as.vector(tmp)

	data$signal = data$V4
	sig = unlist(apply(data[,c("sta","end","signal")], 1, f2, sig="signal", sta="sta", end="end"))
	names(sig) = tmp
	
	#recover the experiment files
	information = strsplit(names(sig) ,"_")
	chr_name = lapply(information, function(x) as.vector(x[1]))
	
	xx = lapply(information, function(x) as.vector(x[2]))
	new_start = lapply(information, function(x) as.vector(as.numeric(x[2])*100 +1))
	
	
	chr_name = as.vector(unlist(chr_name))
	new_start = as.vector(unlist(new_start))
	new_end = new_start + 99
	signal = as.vector(sig)
	
	histone = cbind(chr_name, new_start, new_end,signal)
	histone = as.data.frame(histone)
	colnames(histone) = c("V1","V2","V3","signal")
	histone$ID = paste0(histone$V1,"_",histone$V2,"_",histone$V3)
	
	positive$ID = paste0(positive$V1,"_",positive$V2,"_",positive$V3)

	row.names(histone) = histone$ID
	row.names(positive) = positive$ID
	
	com = intersect(row.names(histone),row.names(positive))
	
	histone = histone[com,]
	
	myoutf = paste0(mytardir,my_nam[i],"_ENCFF172DEA_positive.Rda")
	save(histone,file = myoutf)
}

[4]begin to test the negative counts
rm(list=ls())
#define the function
f1 <- function(x, chr, sta, end) {paste(x[chr], x[sta]:x[end], sep="_")} #separate each base pair
f2 <- function(x, sig, sta, end) 
	{
		time = x[end]-x[sta]+1
		rep(as.character(x[sig]), time)
	}
	
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_counts/"
files = list.files(mydir)

myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_negative.Rda"
load(myinf)
bin.size = 100

mytardir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_counts_bins/"
my_nam = gsub(".bedgraph","",files)

for(i in 2 : length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir,files[i])
	data = read.table(tmpinf,sep="\t",quote=NULL)
	data$V2 = data$V2 +1 
	data$V3 = data$V3 +1 
	data$width =data$V3 - data$V2
	tag = data$width >= 100
	data = data[tag,]
	
	tar = paste0("chr",c(seq(1,22,1),"X","Y"))
	tag = which(data$V1 %in% tar)
	data = data[tag,]
	
	data$sta = ceiling((data$V2)/bin.size)
	data$end = floor((data$V3)/bin.size)
	
	tmp = unlist(apply(data, 1, f1, chr="V1", sta="sta", end="end"))
	tmp = as.vector(tmp)

	data$signal = data$V4
	sig = unlist(apply(data[,c("sta","end","signal")], 1, f2, sig="signal", sta="sta", end="end"))
	names(sig) = tmp
	
	#recover the experiment files
	information = strsplit(names(sig) ,"_")
	chr_name = lapply(information, function(x) as.vector(x[1]))
	
	xx = lapply(information, function(x) as.vector(x[2]))
	new_start = lapply(information, function(x) as.vector(as.numeric(x[2])*100 +1))
	
	
	chr_name = as.vector(unlist(chr_name))
	new_start = as.vector(unlist(new_start))
	new_end = new_start + 99
	signal = as.vector(sig)
	
	histone = cbind(chr_name, new_start, new_end,signal)
	histone = as.data.frame(histone)
	colnames(histone) = c("V1","V2","V3","signal")
	histone$ID = paste0(histone$V1,"_",histone$V2,"_",histone$V3)
	
	negative$ID = paste0(negative$V1,"_",negative$V2,"_",negative$V3)

	row.names(histone) = histone$ID
	row.names(negative) = negative$ID
	
	com = intersect(row.names(histone),row.names(negative))
	
	histone = histone[com,]
	
	myoutf = paste0(mytardir,my_nam[i],"_ENCFF172DEA_negative.Rda")
	save(histone,file = myoutf)
}

[4]begin to model
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_counts_bins/"
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/Histone_file_annotation.txt"

files = list.files(mydir)
name = gsub(".Rda","",files)
name = gsub("_ENCFF172DEA_positive","",name)
name = gsub("_ENCFF172DEA_negative","",name)
name = unique(name)
info = read.table(myinf1,sep="\t",quote=NULL,stringsAsFactors=F)
row.names(info) = info$file_accession
info = info[name,]
info = info[,c("target","file_accession")]

library(randomForest)
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Multi_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Uni_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/moveme.R")

AUC_score = list()

for(i in 1 : nrow(info))
{
	cat("\r",i)
	
	tmpinf1 = paste0(mydir,row.names(info)[i],"_ENCFF172DEA_negative.Rda")
	load(tmpinf1)
	negative = histone
	tag = sample(seq(1,nrow(negative)),10000)
	negative = negative[tag,]
	
	tmpinf2 = paste0(mydir,row.names(info)[i],"_ENCFF172DEA_positive.Rda")
	load(tmpinf2)
	positive = histone
	tag = sample(seq(1,nrow(positive)),10000)
	positive = positive[tag,]
	
	negative$tag = rep(0,nrow(negative))
	positive$tag = rep(1,nrow(positive))
	
	data = rbind(positive,negative)
	data = data[,c("tag","signal")]
	data$signal = as.numeric(as.vector(data$signal))
	data$signal = sign(data$signal)
	data$tag = factor(data$tag, levels =c(0,1))
	fit = ROC_uni(data)
	
	if(mean(fit[[3]])>0.5)
	{
		AUC_score[[i]] = fit[[3]]
	}else{
	
		AUC_score[[i]] = 1-fit[[3]]
	}
}
names(AUC_score) = info$target
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_Histone_AUC.Rda"
save(AUC_score,file = myoutf)

sapply(AUC_score,function(x) mean(x))

 H3K4me3   H3K4me2   H3K4me3   H3K4me1  H3K27me3  H3K27me3  H3K79me2   H3K4me3 
0.6016000 0.6369500 0.6179500 0.5513500 0.5080128 0.5126000 0.5529500 0.5562500 
  H3K9me3   H3K4me3    H3K9ac  H3K36me3   H3K4me3  H3K36me3   H3K9me3  H3K79me2 
0.5077539 0.5375000 0.6042000 0.5229000 0.5315000 0.5293000 0.5228000 0.5444500 
   H3K9ac   H3K9me3     H2AFZ  H3K27me3   H3K4me2   H3K27ac   H3K4me1   H3K4me3 
0.6195000 0.5148000 0.6326500 0.5104289 0.6394500 0.6141500 0.5800500 0.6137000 
    H2AFZ  H4K20me1  H3K36me3  H3K27me3  H4K20me1   H3K27ac  H3K36me3 
0.6481000 0.5056836 0.5246000 0.5086962 0.5148500 0.6134500 0.5212000 

[5]Being to model as all

rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_counts_bins/"
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/Histone_file_annotation.txt"

files = list.files(mydir)
name = gsub(".Rda","",files)
name = gsub("_ENCFF172DEA_positive","",name)
name = gsub("_ENCFF172DEA_negative","",name)
name = unique(name)
info = read.table(myinf1,sep="\t",quote=NULL,stringsAsFactors=F)
row.names(info) = info$file_accession
info = info[name,]
info = info[,c("target","file_accession")]

library(randomForest)
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Multi_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Uni_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/moveme.R")

AUC_score = list()

tag_neg = sample(seq(1,803731),5000)
tag_pos = sample(seq(1,803731),5000)

res = matrix(0, 10000, nrow(info))
colnames(res) = row.names(info)
res = as.data.frame(res)

for(i in 1 : nrow(info))
{
	cat("\r",i)
	
	tmpinf1 = paste0(mydir,row.names(info)[i],"_ENCFF172DEA_negative.Rda")
	load(tmpinf1)
	negative = histone
	negative = negative[tag_neg,]
	
	tmpinf2 = paste0(mydir,row.names(info)[i],"_ENCFF172DEA_positive.Rda")
	load(tmpinf2)
	positive = histone
	positive = positive[tag_pos,]
	
	negative$tag = rep(0,nrow(negative))
	positive$tag = rep(1,nrow(positive))
	
	data = rbind(positive,negative)
	
	if(i == 1)
	{
		row.names(res) = row.names(data)
		res[,i] = data$signal
	}else{
	
		res[,i] = data$signal
	}

}
raw.res = res
res = apply(res,2,function(x) as.numeric(as.vector(x)))
res = sign(res)
tag = c(rep(1,5000),rep(0,5000))
data = cbind(tag,res)
data = as.data.frame(data)
data$tag = factor(data$tag, levels =c(0,1))
fit = ROC(data)
	
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_Histone_integration_AUC.Rda"
save(fit,file = myoutf)
