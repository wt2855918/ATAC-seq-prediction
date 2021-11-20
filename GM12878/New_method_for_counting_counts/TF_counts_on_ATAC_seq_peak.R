[1]Convert the bigwig files to bedgraph
rm(list=ls())
work_mydir  = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Convert_TF_from_bigWig_to_bed/"
dir.create(work_mydir)
setwd(work_mydir)

sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
sam = gsub(".bigWig","",bam_files)

myoutdir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_bed/"
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
	
	curLine = c("cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/")
	writeLines(curLine, conOut)
	
	input_file = bam_files[k]
	
	output = paste0(myoutdir,sam[k],".bedGraph")
	
	curLine = paste0("/lorax/chenglab/yanding/Software/bigWigToBedGraph ", input_file," ", output)
	
	writeLines(curLine, conOut)	
	close(conOut)
}

#begin to submit job
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Convert_TF_from_bigWig_to_bed/"
files = list.files(myDir1)

for(i in 1 : length(files))
{
	cat("\r",i)
	myinf = paste0(myDir1,files[i])
	command = paste0("qsub ", myinf)
	system(command)
}

[2]Begin to calculate the signal for each bin
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/counting_TF_genome/"
dir.create(myDir1)
mytemplate = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Template/Template_for_counting_TF_coverage_genome.R"
mysub = paste(myDir1, "submit.sp", sep="")
setwd(myDir1)

conIn = file(mytemplate, "r")
rawdata = readLines(conIn)
close(conIn)

fnum = 469
for(k in 1:fnum)
{
	data = rawdata
	se = grep("myk = kkkk", data)
	tmp = data[se]
	tmp = gsub("kkkk", k, tmp)
	data[se] = tmp
	myoutf1 = paste("job", k, ".sp", sep="")
	conOut = file(myoutf1, "w")
	writeLines(data, conOut)
	close(conOut)
}



#begin to submit job
work_mydir  = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/counting_TF_genome/"
setwd(work_mydir)
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/counting_TF_genome/"
files = list.files(myDir1)

target = paste0("job",seq(1,20,1),".sp")

for(i in 1 : length(files))
{
	cat("\r",i)
	myinf = paste0(myDir1,files[i])
	command = paste0("qsub ", myinf)
	system(command)
}


[3]Begin to sort the positive and negative files
rm(list=ls())
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_positive.Rda"
load(myinf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/ATAC_seq_peak_processed_bed/ENCFF172DEA_positive.txt"
write.table(positive,myoutf,sep="\t",quote=F,col.names=F,row.names=F)

rm(list=ls())
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_negative.Rda"
load(myinf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/ATAC_seq_peak_processed_bed/ENCFF172DEA_negative.txt"
write.table(negative,myoutf,sep="\t",quote=F,col.names=F,row.names=F)

cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/ATAC_seq_peak_processed_bed/
cut -f 1,2,3 ENCFF172DEA_positive.txt > ENCFF172DEA_positive.bed
cut -f 1,2,3 ENCFF172DEA_negative.txt > ENCFF172DEA_negative.bed

/lorax/chenglab/yanding/Software/bedtools2/bin/sortBed -i ENCFF172DEA_positive.bed > ENCFF172DEA_positive_sort.bed
/lorax/chenglab/yanding/Software/bedtools2/bin/sortBed -i ENCFF172DEA_negative.bed > ENCFF172DEA_negative_sort.bed

[4]Calculate the positive counts
rm(list=ls())
work_mydir  = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF_pos/"
dir.create(work_mydir)
setwd(work_mydir)

sourcedir = "/lorax/chenglab/yanding/m6A_integration/Data/HepG2/TF_bed_sort/"
bam_files = list.files(sourcedir)
sam = gsub("_sort.bedGraph","",bam_files)

myoutdir = "/lorax/chenglab/yanding/m6A_integration/Data/HepG2/RBP_logFC_m6A_pos_strand/"
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
	
	curLine = c("cd /lorax/chenglab/yanding/m6A_integration/Data/HepG2/RBP_logFC_sort_bed_pos/")
	writeLines(curLine, conOut)
	
	input_file2 = bam_files[k]
	input_file1 = "/lorax/chenglab/yanding/m6A_integration/Raw_data/HepG2/m6A/HepG2_pos_sort_pos_strand.bed"
	
	output = paste0(myoutdir,sam[k],"_pos.bedGraph")
	
	curLine = paste0("/lorax/chenglab/yanding/Software/bedtools2/bin/mapBed -a ", input_file1," -b ",input_file2," -c 4 -o mean > ",output)
	
	writeLines(curLine, conOut)	
	close(conOut)
}

#begin to submit job
myDir1 = "/lorax/chenglab/yanding/m6A_integration/Scripts/HepG2/Counting_RBP_pos_pos_strand/"
files = list.files(myDir1)

for(i in 1 : length(files))
{
	cat("\r",i)
	myinf = paste0(myDir1,files[i])
	command = paste0("qsub ", myinf)
	system(command)
}




#PBS -l walltime=200:00:00 -l feature=celln
/usr/bin/R --vanilla <<EOF

rm(list=ls())
#define the function
f1 <- function(x, chr, sta, end) {paste(x[chr], x[sta]:x[end], sep="_")} #separate each base pair
f2 <- function(x, sig, sta, end) 
	{
		time = x[end]-x[sta]+1
		rep(as.character(x[sig]), time)
	}
	
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]


myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_negative.Rda"
load(myinf)
bin.size = 100

mytardir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_counts_bins_logFC/"
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
	
myoutf = paste0(mytardir,my_nam[i],"_ENCFF172DEA_negative.Rda")
save(histone,file = myoutf)
EOF




