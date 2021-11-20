[1]Data download
cd /lorax/chenglab/yanding/Pub_Dat/GRch38/Fastq_files/
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.8.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.9.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.16.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.18.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz 

[2]begin to process ATAC-seq peaks batch 1
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/ATAC_seq/Batch_1/"
myDir2 = "/lorax/chenglab/yanding/Pub_Dat/GRch38/Fastq_files/"

mytarDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_sequence/Batch_1/"
dir.create(mytarDir1)

mychr = paste("chr", c(1:22, "X"),  sep="")
files = list.files(myDir1)
tag = grep(".bed",files)
files = files[tag]

mytag = gsub(".bed", "", files)
#myk = 3#kkkk
myinf1 = paste(myDir1, files, sep="")

for(k in 1:length(mytag))
{
	tmpf = paste(mytarDir1, mytag[k], "_peak.fa", sep="")
	conOut = file(tmpf, "w")
	close(conOut)
}


for(k in 1 : length(mytag))
{
	
	myp = read.table(myinf1[k], sep="\t", stringsAsFactors=F)
	#myp = myp[order(myp$V9,decreasing=T),]
	#myp = myp[1:500,]
	#id = paste0(myp$V1,"_",myp$V2,"_",myp$V3)
	#myp$id = id
	myp = myp[,c(1,2,3)]
	myp = unique(myp)

	colnames(myp)[1:3]= c("chr","start","end")
	#myp[,"length"] = myp[,"end"]-myp[,"start"]
	#myp[,"abs_summit"] = myp[,"start"] + floor(myp[,"length"]/2)

	#mylen = myp[, "length"]
	#xx = myp[, "abs_summit"]-100
	#mysta = ifelse(xx>myp[, "start"], xx, myp[, "start"])
	#xx = myp[, "abs_summit"]+100
	#myend = ifelse(xx<myp[, "end"], xx, myp[, "end"])
	#tag = mylen>=501
	#myp[tag==1,2] = mysta[tag==1]
	#myp[tag==1,3] = myend[tag==1]
	#myp = myp[mylen>20,]
	tmpf = paste(mytarDir1, mytag[k], "_peak.fa", sep="")
	conOut = file(tmpf, "w")

	summary = list()
	bad_id = list()

	for(p in 1:length(mychr))
	{
		cat("\r",k,"-->,",p)
		mytmpf = paste(myDir2,  mychr[p], ".fa", sep="")
		conIn = file(mytmpf, "r")
		myseq= readLines(conIn, -1)
		close(conIn)
		myseq =myseq[-1]
		myseq = paste(myseq, collapse="")  
	
		tag = myp[,1]==mychr[p]
		if(sum(tag)==0) next
		
		curp = myp[tag==1,]
		
		chrom_stat = rep(0,nrow(curp))
		ID = rep(0,nrow(curp))
		
			for(i in 1:nrow(curp))
		  	{
			  tmpseq = toupper(substr(myseq, curp[i,2], curp[i,3]))
			  
			  test = strsplit(tmpseq,"")[[1]]
			  
			  	tag = grep("N",test)
			  	if(length(tag)>0)
			  	{
			  		chrom_stat[i] = 1
			  		ID[i] = paste0(curp$chr[i],"_",curp$start[i],"_",curp$end[i])
			  	}else{
			  		writeLines(paste(">",curp[i,1], "_", curp[i,2], "_", curp[i,3], sep=""), conOut)
			  		writeLines(tmpseq, conOut)
			  	}  
		 	 }
	
		summary[[p]] = chrom_stat
		bad_id[[p]] = ID
	}
	close(conOut)  	

	#seq_test = unlist(summary)#only 7 peaks have N
	#peak_to_delete = unique(unlist(bad_id))
	#peak_to_delete = peak_to_delete[2:length(peak_to_delete)]

	#k = myk =1
	#myp = read.table(myinf1[k], sep="\t", stringsAsFactors=F)
	#id = paste0(myp$V1,"_",myp$V2,"_",myp$V3)
	#tag = which(id %in% peak_to_delete)
	#myp = myp[-tag,]

	#myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_processed_peaks/GM12878_batch_1_peaks.txt"
	#write.table(myp,myoutf,sep="\t",quote=F)
}
[2]Begin to check
rm(list=ls())
library("Biostrings")
library("reticulate")
source_python("/lorax/chenglab/yanding/ATAC_seq_integration/Function/KMER_Counts.py")

myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_sequence/Batch_1/"
files = list.files(myinf)
my_nam = gsub(".fa","",files)

for(k in 1 : length(files))
{
	tmpinf = paste0(myinf, files[k])
	fastaFile <- readDNAStringSet(tmpinf)
	seq_name = names(fastaFile)
	sequence = paste(fastaFile)
	df <- data.frame(seq_name, sequence)

	#numFeature = 4096 (4^6)
	for(i in 1 : length(sequence))
	{
		cat("\r",k,"-->",i)	
		features = count_kmers(sequence[i], 4096)
	
		if(i == 1)
		{
			res = as.data.frame(features)
		}else{
			res = cbind(res,features)
		}	
}

	row.names(res) = paste0("Sig",seq(1,4096,1))
	colnames(res) = seq_name

	myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_sequence_feature/Batch_1/",my_nam[k],".Rda")
	save(res,file = myoutf)
}

[4]Begint o model self prediction and 
library(neuralnet)
#myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_processed_peaks/GM12878_batch_1_peaks.txt"
mydir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/ATAC_seq/Batch_1/"
mydir2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_sequence_feature/Batch_1/"

files = list.files(mydir1)
files = files[c(2:length(files))] #remove the first combined file
my_nam = gsub(".bed","",files)

for(i in 1 : length(files))
{
	cat("\r",i)
	
	training_nam = my_nam[i]
	
	tmpinf1 = paste0(mydir1,training_nam,".bed")
	tmpinf2 = paste0(mydir2,training_nam,"_peak.Rda")
	
	training_data = read.table(tmpinf1,sep="\t",quote=NULL)
	load(tmpinf2)
	
	ID = paste0(training_data$V1,"_",training_data$V2,"_",training_data$V3)
	training_data$ID = ID
	training_data = aggregate(training_data[, 7:9], list(training_data$ID), mean)
	colnames(training_data)[1] = "ID"
	row.names(training_data) = training_data$ID

	com = intersect(row.names(training_data), colnames(res))
	training_data = training_data[com,]
	res = res[,com]
	res = t(res)
	
	#training_data$Value = (training_data$V7 - mean(training_data$V7))/sd(training_data$V7)
	training_data$Value  = training_data$V7
	df = cbind(Value = training_data$Value, res)
	
	nn=neuralnet(Value~.,data=df, hidden=100, act.fct = "logistic",
                linear.output = T, algorithm = "backprop",learningrate =0.0001)
    
    self_validation = compute(nn, res)
    predicted_value = self_validation$net.result
    predicted_value = (predicted_value - mean(predicted_value))/sd(predicted_value)
    
    #cor(nn$response[,"Value"], df[,"Value"])
    cor(predicted_value, df[,"Value"])
}


softplus <- function(x) log(1+exp(x))
relu <- function(x) sapply(x, function(z) max(0,z))
net$act.fct(x)


myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/ATAC_seq/Batch_1/ENCFF287FYP.bed"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/ATAC_seq/Batch_1/ENCFF348YMZ.bed"
myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/ATAC_seq/Batch_1/ENCFF912ESR.bed"


myinf4 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_sequence_feature/Batch_1_ENCFF172DEA_summary.Rda"

data1 = read.table(myinf1,sep="\t",quote=NULL)
data2 = read.table(myinf2,sep="\t",quote=NULL)
data3 = read.table(myinf3,sep="\t",quote=NULL)

ID1 = paste0(data1$V1,"_",data1$V2,"_",data1$V3)
ID2 = paste0(data2$V1,"_",data2$V2,"_",data2$V3)
ID3 = paste0(data3$V1,"_",data3$V2,"_",data3$V3)

data1$ID = ID1
data2$ID = ID2
data3$ID = ID3

length(intersect(ID1,colnames(res)))
length(intersect(ID2,colnames(res)))
length(intersect(ID3,colnames(res)))

load(myinf2)

ID = paste0(data$V1,"_",data$V2,"_",data$V3)
data$ID = ID
#did not include Y

tag = duplicated(ID)
ID[tag][1]

[5]Begin to create the control region
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/ATAC_seq/Batch_1/"
myDir2 = "/lorax/chenglab/yanding/Pub_Dat/GRch38/Fastq_files/"

mytarDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_sequence_ctrl/Batch_1/"
dir.create(mytarDir1)

files = list.files(myDir1)
my_nam = gsub(".bed","",files)

####################################################
#define the left name function
f1 <- function(x, chr, sta, end, length, upstream) 
{
	new_start = as.numeric(x[sta]) - as.numeric(x[length]) - upstream
	new_end = as.numeric(x[end]) - as.numeric(x[length]) - upstream
	paste(x[chr], new_start, new_end, "L", sep="_")
} 

for(i in 1: length(files))
{
	cat("\r",i)
	tmpinf = paste0(myDir1, files[i])
	data = read.table(tmpinf,sep="\t",quote=NULL,stringsAsFactors=F)
	data = as.data.frame(data)
	data$V2 = as.numeric(data$V2)
	data$V3 = as.numeric(data$V3)
	data$width = data$V3 - data$V2
	data = as.data.frame(data)
	data = data[,c("V1","V2","V3","width")]
	data = unique(data)
	data$ID = paste0(data$V1,"_",data$V2,"_",data$V3)
	left_name = apply(data, 1, f1, chr="V1", sta="V2", end="V3", length = "width", upstream = 1)
	
	information = strsplit(left_name ,"_")
	chr_name = lapply(information, function(x) as.vector(x[1]))
	new_start = lapply(information, function(x) as.vector(x[2]))
	new_end = lapply(information, function(x) as.vector(x[3]))
	
	chr_name = as.vector(unlist(chr_name))
	new_start = as.vector(unlist(new_start))
	new_end = as.vector(unlist(new_end))
	
	left_control = cbind(chr_name, new_start, new_end)
	left_control = as.data.frame(left_control)
	
	myoutf = paste0(mytarDir1,my_nam[i],"_left.txt")
	write.table(left_control,myoutf,sep="\t",quote=F)
}
########################################################################################################
#define the right name function
f1 <- function(x, chr, sta, end, length, dnstream) 
{
	new_start = as.numeric(x[sta]) + as.numeric(x[length]) + dnstream
	new_end = as.numeric(x[end]) + as.numeric(x[length]) + dnstream
	paste(x[chr], new_start, new_end, "R", sep="_")
} 

for(i in 1: length(files))
{
	cat("\r",i)
	tmpinf = paste0(myDir1, files[i])
	data = read.table(tmpinf,sep="\t",quote=NULL,stringsAsFactors=F)
	data = as.data.frame(data)
	data$V2 = as.numeric(data$V2)
	data$V3 = as.numeric(data$V3)
	data$width = data$V3 - data$V2
	data = as.data.frame(data)
	data = data[,c("V1","V2","V3","width")]
	data = unique(data)
	data$ID = paste0(data$V1,"_",data$V2,"_",data$V3)
	right_name = apply(data, 1, f1, chr="V1", sta="V2", end="V3", length = "width", dnstream = 1)
	
	information = strsplit(right_name ,"_")
	chr_name = lapply(information, function(x) as.vector(x[1]))
	new_start = lapply(information, function(x) as.vector(x[2]))
	new_end = lapply(information, function(x) as.vector(x[3]))
	
	chr_name = as.vector(unlist(chr_name))
	new_start = as.vector(unlist(new_start))
	new_end = as.vector(unlist(new_end))
	
	right_control = cbind(chr_name, new_start, new_end)
	right_control = as.data.frame(right_control)
	
	myoutf = paste0(mytarDir1,my_nam[i],"_right.txt")
	write.table(right_control,myoutf,sep="\t",quote=F)
}

[6]Begin to examine the overlap
rm(list=ls())
mydir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/ATAC_seq/Batch_1/"
mydir2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_sequence_ctrl/Batch_1/"

files = list.files(mydir1)

my_nam = gsub(".bed","",files)

#define function
f1 <- function(x, chr, sta, end) {paste(x[chr], x[sta]:x[end], sep="_")} #separate each base pair
f2 <- function(x, sig, sta, end) 
	{
		time = x[end]-x[sta]+1
		rep(as.character(x[sig]), time)
	}
bin.size = 100

for(i in 1 : length(my_nam))
{
	cat("\r",i)
	
	tmpinf1 = paste0(mydir1,my_nam[i],".bed")
	tmpinf2 = paste0(mydir2,my_nam[i],"_left.txt")
	tmpinf3 = paste0(mydir2,my_nam[i],"_right.txt")
	
	res1 = read.table(tmpinf1,sep="\t",quote=NULL,stringsAsFactors=F)
	res2 = read.table(tmpinf2,sep="\t",quote=NULL,stringsAsFactors=F)
	res3 = read.table(tmpinf3,sep="\t",quote=NULL,stringsAsFactors=F)
	
	###################################
	#Process the experimental group
	####################################
	res1 = res1[,c("V1","V2","V3")]
	res1 = unique(res1)
	res1$ID = paste0(res1$V1,"_",res1$V2,"_",res1$V3)
	
	res1$sta = ceiling((res1$V2)/bin.size)
	res1$end = floor((res1$V3)/bin.size)
	tmp1 = unlist(apply(res1, 1, f1, chr="V1", sta="sta", end="end"))
	tmp1 = as.vector(tmp1)

	res1$signal = seq(1,nrow(res1),1)
	#res1$signal = res1$ID
	sig1 = unlist(apply(res1[,c("sta","end","signal")], 1, f2, sig="signal", sta="sta", end="end"))
	names(sig1) = tmp1
	peak_1 = sig1
	
	#recover the experiment files
	information = strsplit(names(peak_1) ,"_")
	chr_name = lapply(information, function(x) as.vector(x[1]))
	new_start = lapply(information, function(x) as.vector(as.numeric(x[2])*100 +1))
	
	
	chr_name = as.vector(unlist(chr_name))
	new_start = as.vector(unlist(new_start))
	new_end = new_start + 99
	
	positive = cbind(chr_name, new_start, new_end)
	positive = as.data.frame(positive)
	colnames(positive) = c("V1","V2","V3")
	positive = unique(positive)
	
	myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/",my_nam[i],"_positive.Rda")
	save(positive, file =myoutf)
	
	###################################
	#Process the left control
	####################################
	res2 = res2[,c("chr_name","new_start","new_end")]
	res2 = unique(res2)
	res2$ID = paste0(res2$chr_name,"_",res2$new_start,"_",res2$new_end)
	
	res2$sta = ceiling((res2$new_start)/bin.size)
	res2$end = floor((res2$new_end)/bin.size)
	tmp2 = unlist(apply(res2, 1, f1, chr="chr_name", sta="sta", end="end"))
	tmp2 = as.vector(tmp2)

	res2$signal = seq(1,nrow(res2),1)
	#res2$signal = res2$ID
	sig2 = unlist(apply(res2[,c("sta","end","signal")], 1, f2, sig="signal", sta="sta", end="end"))
	names(sig2) = tmp2
	peak_2 = sig2
	
	###################################
	#Process the right control
	####################################
	res3 = res3[,c("chr_name","new_start","new_end")]
	res3 = unique(res3)
	res3$ID = paste0(res3$chr_name,"_",res3$new_start,"_",res3$new_end)
	
	res3$sta = ceiling((res3$new_start)/bin.size)
	res3$end = floor((res3$new_end)/bin.size)
	tmp3 = unlist(apply(res3, 1, f1, chr="chr_name", sta="sta", end="end"))
	tmp3 = as.vector(tmp3)

	res3$signal = seq(1,nrow(res3),1)
	#res2$signal = res2$ID
	sig3 = unlist(apply(res3[,c("sta","end","signal")], 1, f2, sig="signal", sta="sta", end="end"))
	names(sig3) = tmp3
	peak_3 = sig3
	
	#begin to examine the overlap
	
	com = intersect(names(peak_1),names(peak_2))
	tag = which(names(peak_2) %in% com)
	left_peak_kept = peak_2[-tag]
	
	com = intersect(names(peak_1),names(peak_3))
	tag = which(names(peak_3) %in% com)
	right_peak_kept = peak_3[-tag]

	peak_kept = c(left_peak_kept, right_peak_kept)
	peak_kept = names(peak_kept)
	peak_kept = unique(peak_kept)
	
	#chose the exact same number
	tag = sample(1:length(peak_kept), nrow(positive), replace=F)
	peak_kept = peak_kept[tag]
	
	#recover the experiment files
	control = strsplit(peak_kept ,"_")
	chr_name = lapply(control, function(x) as.vector(x[1]))
	new_start = lapply(control, function(x) as.vector(as.numeric(x[2])*100 +1))
	
	chr_name = as.vector(unlist(chr_name))
	new_start = as.vector(unlist(new_start))
	new_end = new_start + 99
	
	negative = cbind(chr_name, new_start, new_end)
	negative = as.data.frame(negative)
	colnames(negative) = c("V1","V2","V3")
	negative = unique(negative)
	
	myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/",my_nam[i],"_negative.Rda")
	save(negative, file =myoutf)
	
}
[7]extract the sequence again #template
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/"
myDir2 = "/lorax/chenglab/yanding/Pub_Dat/GRch38/Fastq_files/"

mytarDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1/"
dir.create(mytarDir1)

mychr = paste("chr", c(1:22, "X"),  sep="")
files = list.files(myDir1)
tag = grep(".Rda",files)
files = files[tag]

mytag = gsub(".Rda", "", files)
#myk = 3#kkkk
myinf1 = paste(myDir1, files, sep="")

for(k in 1:length(mytag))
{
	tmpf = paste(mytarDir1, mytag[k], "_peak.fa", sep="")
	conOut = file(tmpf, "w")
	close(conOut)
}


for(k in 1 : length(mytag))
{
	
	x = load(myinf1[k])
	myp = get(x)
	myp = as.data.frame(myp,stringsAsFactors=F)
	rm(x)

	colnames(myp)[1:3]= c("chr","start","end")
	myp[,"chr"] = as.character(as.vector(myp[,"chr"]))
	myp[,"start"] = as.numeric(as.vector(myp[,"start"]))
	myp[,"end"] = as.numeric(as.vector(myp[,"end"]))
	
	tmpf = paste(mytarDir1, mytag[k], "_peak.fa", sep="")
	conOut = file(tmpf, "w")

	summary = list()
	bad_id = list()

	for(p in 1:length(mychr))
	{
		
		mytmpf = paste(myDir2,  mychr[p], ".fa", sep="")
		conIn = file(mytmpf, "r")
		myseq= readLines(conIn, -1)
		close(conIn)
		myseq =myseq[-1]
		myseq = paste(myseq, collapse="")  
	
		tag = myp[,1]==mychr[p]
		if(sum(tag)==0) next
		
		curp = myp[tag==1,]
		
		chrom_stat = rep(0,nrow(curp))
		ID = rep(0,nrow(curp))
		
			for(i in 1:nrow(curp))
		  	{
		  		cat("\r",k,"-->",p,"-->",i)	
			  tmpseq = toupper(substr(myseq, curp[i,2], curp[i,3]))
			  
			  test = strsplit(tmpseq,"")[[1]]
			  
			  	tag = grep("N",test)
			  	if(length(tag)>0)
			  	{
			  		chrom_stat[i] = 1
			  		ID[i] = paste0(curp$chr[i],"_",curp$start[i],"_",curp$end[i])
			  	}else{
			  		writeLines(paste(">",curp[i,1], "_", curp[i,2], "_", curp[i,3], sep=""), conOut)
			  		writeLines(tmpseq, conOut)
			  	}  
		 	 }
	
		#summary[[p]] = chrom_stat
		#bad_id[[p]] = ID
	}
	close(conOut)  	
}

[8]Use template to get the sequence
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Sequence_get/"
mytemplate = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Template/Template_for_getting_sequence.txt"
mysub = paste(myDir1, "submit.sp", sep="")
setwd(myDir1)

conIn = file(mytemplate, "r")
rawdata = readLines(conIn)
close(conIn)

fnum = 8
for(k in 1:fnum)
{
	data = rawdata
	se = grep("k = mykkk", data)
	tmp = data[se]
	tmp = gsub("mykkk", k, tmp)
	data[se] = tmp
	myoutf1 = paste("job", k, ".sp", sep="")
	conOut = file(myoutf1, "w")
	writeLines(data, conOut)
	close(conOut)
}


[9]Count Kmers template
rm(list=ls())
library("Biostrings")
library("reticulate")
source_python("/lorax/chenglab/yanding/ATAC_seq_integration/Function/KMER_Counts.py")

myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1/"
files = list.files(myinf)
my_nam = gsub(".fa","",files)

for(k in 1 : length(files))
{
	tmpinf = paste0(myinf, files[k])
	fastaFile <- readDNAStringSet(tmpinf)
	seq_name = names(fastaFile)
	sequence = paste(fastaFile)
	#df <- data.frame(seq_name, sequence)

	#numFeature = 4096 (4^6)
	for(i in 1 : length(sequence))
	{
		cat("\r",k,"-->",i)	
		features = count_kmers(sequence[i], 4096)
	
		if(i == 1)
		{
			res = as.data.frame(features)
		}else{
			res = cbind(res,features)
		}	
	}

	row.names(res) = paste0("Sig",seq(1,4096,1))
	colnames(res) = seq_name

	myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence_feature/Batch_1/",my_nam[k],".Rda")
	save(res,file = myoutf)
}

[10]Count kmers job
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Kmer_get/"
dir.create(myDir1)
mytemplate = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Template/Template_for_counting_kmers.txt"
mysub = paste(myDir1, "submit.sp", sep="")
setwd(myDir1)

conIn = file(mytemplate, "r")
rawdata = readLines(conIn)
close(conIn)

fnum = 8
for(k in 1:fnum)
{
	data = rawdata
	se = grep("k = mykkk", data)
	tmp = data[se]
	tmp = gsub("mykkk", k, tmp)
	data[se] = tmp
	myoutf1 = paste("job", k, ".sp", sep="")
	conOut = file(myoutf1, "w")
	writeLines(data, conOut)
	close(conOut)
}



[9]Random forest

rm(list=ls())
library(randomForest)
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence_feature/Batch_1/ENCFF172DEA_positive_peak.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence_feature/Batch_1/ENCFF172DEA_negative_peak.Rda"
load(myinf1)

label = sample(seq(1,803717,1),10000)

pos_sam = res[,label]
write.table(pos_sam, "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence_feature/Batch_1/ENCFF172DEA_trans_sample_positive_peak.txt", sep="\t",quote=F)

rm(res)
load(myinf2)
#res = unique(res)


label = sample(seq(1,803691,1),10000)

neg_sam = res[,label]
write.table(neg_sam, "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence_feature/Batch_1/ENCFF172DEA_trans_sample_negative_peak.txt", sep="\t",quote=F)



rm(list=ls())
library(randomForest)
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence_feature/Batch_1/ENCFF172DEA_positive_peak.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence_feature/Batch_1/ENCFF172DEA_negative_peak.Rda"
load(myinf1)

label = sample(seq(1,803717,1),20000)

pos_sam = res[,label]
write.table(pos_sam, "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence_feature/Batch_1/ENCFF172DEA_trans_sample_positive_peak_20k.txt", sep="\t",quote=F)

rm(res)
load(myinf2)
#res = unique(res)


label = sample(seq(1,803691,1),20000)

neg_sam = res[,label]
write.table(neg_sam, "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence_feature/Batch_1/ENCFF172DEA_trans_sample_negative_peak_20k.txt", sep="\t",quote=F)






