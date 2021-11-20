[1]Calcualte the coverage
rm(list=ls())
work_mydir  = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF/"
dir.create(work_mydir)
setwd(work_mydir)

sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq/"
bam_files = list.files(sourcedir)
tag = grep("bai",bam_files)
bam_files = bam_files[-tag]
sam = gsub(".bam","",bam_files)

myoutdir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_counts/"
dir.create(myoutdir)

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
	
	curLine = c("cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq/")
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
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF/"
files = list.files(myDir1)

for(i in 1 : length(files))
{
	cat("\r",i)
	myinf = paste0(myDir1,files[i])
	command = paste0("qsub ", myinf)
	system(command)
}

[2]Begin to test the positive counts
rm(list=ls())
#define the function
f1 <- function(x, chr, sta, end) {paste(x[chr], x[sta]:x[end], sep="_")} #separate each base pair
f2 <- function(x, sig, sta, end) 
	{
		time = x[end]-x[sta]+1
		rep(as.character(x[sig]), time)
	}
	
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_counts/"
files = list.files(mydir)

myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_positive.Rda"
load(myinf)
bin.size = 100

mytardir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_counts_bins/"
dir.create(mytardir)
my_nam = gsub(".bedgraph","",files)

for(i in 1 : length(files))
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
	
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_counts/"
files = list.files(mydir)

myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_negative.Rda"
load(myinf)
bin.size = 100

mytardir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_counts_bins/"
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