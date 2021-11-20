/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/Prostate/ATAC_seq
[1]Process peak files
rm(list=ls())
library(liftOver)
library(tidyverse)
library(GenomicRanges)
library(ChIPpeakAnno)
library(rtracklayer)

myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/Prostate/ATAC_seq/"
myDir2 = "/lorax/chenglab/yanding/Pub_Dat/GRch38/Fastq_files/"

mytarDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/ATAC_seq_peak_sequence_ctrl/"
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
###################################################
#define the convert function
peakDF2GRanges <- function(peak.df) {
    peak.gr=GRanges(seqnames=peak.df[,1],
        ranges=IRanges(peak.df[,2], peak.df[,3]))
    cn <- colnames(peak.df)
    if (length(cn) > 3) {
        for (i in 4:length(cn)) {
            mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
        }
    }
    return(peak.gr)
}

#define the liftover
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)

for(i in 1: length(files))
{
	cat("\r",i)
	tmpinf = paste0(myDir1, files[i])
	data = read.table(tmpinf,sep="\t",quote=NULL,stringsAsFactors=F)
	data = as.data.frame(data)
	data = data[,c("V1","V2","V3")]
	data = unique(data)
	
	#liftover to hg19
	data = peakDF2GRanges(data)
	seqlevelsStyle(data) = "UCSC"
	cur19 = liftOver(data, ch)
	cur19 = unlist(cur19)
	genome(cur19) = "hg19"
	cur19 = new("gwaswloc", cur19)
	
	data = as.data.frame(cur19)
	data = unique(data)
	colnames(data) = c("V1","V2","V3","width")
	
	tag = data$width >100
	data = data[tag,]
	
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
	
	#liftover to hg19
	data = peakDF2GRanges(data)
	seqlevelsStyle(data) = "UCSC"
	cur19 = liftOver(data, ch)
	cur19 = unlist(cur19)
	genome(cur19) = "hg19"
	cur19 = new("gwaswloc", cur19)
	
	data = as.data.frame(cur19)
	data = data[,c("seqnames", "start","end","width")]
	data = unique(data)
	colnames(data) = c("V1","V2","V3","width")
	
	tag = data$width >100
	data = data[tag,]
	
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

[2]Get the positive and negative range
rm(list=ls())
library(liftOver)
library(tidyverse)
library(GenomicRanges)
library(ChIPpeakAnno)
library(rtracklayer)

###################################################
#define the convert function
peakDF2GRanges <- function(peak.df) {
    peak.gr=GRanges(seqnames=peak.df[,1],
        ranges=IRanges(peak.df[,2], peak.df[,3]))
    cn <- colnames(peak.df)
    if (length(cn) > 3) {
        for (i in 4:length(cn)) {
            mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
        }
    }
    return(peak.gr)
}

#define the liftover
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)



mydir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/Prostate/ATAC_seq/"
mydir2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/ATAC_seq_peak_sequence_ctrl/"

files = list.files(mydir1)

my_nam = gsub(".bed","",files)

#define function
f1 <- function(x, chr, sta, end) {paste(x[chr], x[sta]:x[end], sep="_")} #separate each base pair
f2 <- function(x, sig, sta, end) 
	{
		time = x[end]-x[sta]+1
		rep(as.character(x[sig]), time)
	}
bin.size = 200

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
	
	#liftover to hg19
	res1 = peakDF2GRanges(res1)
	seqlevelsStyle(res1) = "UCSC"
	cur19 = liftOver(res1, ch)
	cur19 = unlist(cur19)
	genome(cur19) = "hg19"
	cur19 = new("gwaswloc", cur19)
	res1 = as.data.frame(cur19)
	res1 = unique(res1)
	
	tag = res1$width >200
	res1 = res1[tag,]
	
	res1 = res1[,c("seqnames","start","end")]
	colnames(res1) = c("V1","V2","V3")
	
	res1$ID = paste0(res1$V1,"_",res1$V2,"_",res1$V3)
	
	res1$sta = ceiling((res1$V2+1)/bin.size)
	res1$end = ceiling((res1$V3+1)/bin.size)
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
	new_start = lapply(information, function(x) as.vector(as.numeric(x[2])*100))
	
	
	chr_name = as.vector(unlist(chr_name))
	new_start = as.vector(unlist(new_start))
	new_end = new_start + 199
	
	positive = cbind(chr_name, new_start, new_end)
	positive = as.data.frame(positive)
	colnames(positive) = c("V1","V2","V3")
	positive = unique(positive)
	
	myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/ATAC_seq_peak_processed_bed/",my_nam[i],"_positive.Rda")
	save(positive, file =myoutf)
	
	###################################
	#Process the left control
	####################################
	res2 = res2[,c("chr_name","new_start","new_end")]
	res2 = unique(res2)
	res2$ID = paste0(res2$chr_name,"_",res2$new_start,"_",res2$new_end)
	
	res2$sta = ceiling((res2$new_start+1)/bin.size)
	res2$end = ceiling((res2$new_end+1)/bin.size)
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
	
	res3$sta = ceiling((res3$new_start+1)/bin.size)
	res3$end = ceiling((res3$new_end+1)/bin.size)
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
	new_start = lapply(control, function(x) as.vector(as.numeric(x[2])*100))
	
	chr_name = as.vector(unlist(chr_name))
	new_start = as.vector(unlist(new_start))
	new_end = new_start + 199
	
	negative = cbind(chr_name, new_start, new_end)
	negative = as.data.frame(negative)
	colnames(negative) = c("V1","V2","V3")
	negative = unique(negative)
	
	myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/ATAC_seq_peak_processed_bed/",my_nam[i],"_negative.Rda")
	save(negative, file =myoutf)
	
}

[2]Bin the genome 200bp
rm(list=ls())
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tiles <- tileGenome(seqinfo(txdb), tilewidth=200,
                    cut.last.tile.in.chrom=TRUE)
myoutf = "/lorax/chenglab/yanding/Pub_Dat/Reference_genome/Hg19/Bin_200bp_hg19_genome.Rda"
save(tiles, file = myoutf)

[3]Define the padding region
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/ATAC_seq_peak_processed_bed/ENCFF046YUJ_positive.Rda"
load(myinf1)
tag = sample(1:nrow(positive),20000)
training_pos = positive[tag,]

training_pos$start = as.numeric(as.vector(training_pos$V2)) - 900
training_pos$end = as.numeric(as.vector(training_pos$V3)) + 900
training_pos = training_pos[,c("V1","start","end")]
training_pos$ID = paste0(training_pos$V1,"_",training_pos$start,"_",training_pos$end)

load("/lorax/chenglab/yanding/Pub_Dat/Reference_genome/Hg19/Bin_200bp_hg19_genome.Rda")
negative = as.data.frame(tiles)
negative$start = negative$start - 1
negative$end = negative$end - 1

negative$start = as.numeric(as.vector(negative$start)) - 900
negative$end = as.numeric(as.vector(negative$end)) + 900
negative = negative[,c("seqnames","start","end" )]
negative$ID = paste0(negative$seqnames,"_",negative$start,"_",negative$end)

tag = which(negative$ID %in% training_pos$ID)
negative = negative[-tag,]
tag = negative$start >0 & negative$end >0
negative = negative[tag,]
tag = sample(1:nrow(negative),50000)
negative = negative[tag,]

library(ChIPpeakAnno)

library(GenomicRanges)

peakDF2GRanges <- function(peak.df) {
    peak.gr=GRanges(seqnames=peak.df[,1],
        ranges=IRanges(peak.df[,2], peak.df[,3]))
    cn <- colnames(peak.df)
    if (length(cn) > 3) {
        for (i in 4:length(cn)) {
            mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
        }
    }
    return(peak.gr)
}


res1 = peakDF2GRanges(training_pos)
res2 = peakDF2GRanges(negative)
			
shared_peaks = findOverlapsOfPeaks(res1, res2, connectedPeaks="keepAll")
			
peak = shared_peaks$venn_cnt

peak
     res1 res2 Counts count.res1 count.res2
[1,]    0    0      0          0          0
[2,]    0    1  49607          0      49607
[3,]    1    0  19573      19573          0
[4,]    1    1    390        427        393

all.peaks <- shared_peaks$all.peaks
res1.renamed <- all.peaks$res1
res2.renamed <- all.peaks$res2
res1 = as.data.frame(res1.renamed)
res2 = as.data.frame(res2.renamed)

peakNames <- melt(shared_peaks$peaklist[['res1///res2']]$peakNames, value.name="merged.peak.id")
peakNames <- peakNames$value.value

tag = which(row.names(res1) %in% peakNames)
res1 = res1[-tag,]

tag = which(row.names(res2) %in% peakNames)
res2 = res2[-tag,]

res1 = as.data.frame(res1)
res2 = as.data.frame(res2)

res1 = res1[,c(1,2,3)]
res2 = res2[,c(1,2,3)]
colnames(res1) = colnames(res2) = c("V1","V2","V3")

res1 = peakDF2GRanges(res1)
res2 = peakDF2GRanges(res2)
shared_peaks = findOverlapsOfPeaks(res1, res2, connectedPeaks="keepAll")
			
peak = shared_peaks$venn_cnt

		
training_pos = as.data.frame(res1)
negative = as.data.frame(res2)


myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/ATAC_seq_peak_processed_bed/broad_positive_200bp_ENCFF046YUJ.Rda"
save(training_pos,file =myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/ATAC_seq_peak_processed_bed/broad_negative_200bp_ENCFF046YUJ.Rda"
save(negative,file =myoutf)

[3]Get the sequence
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/ATAC_seq_peak_processed_bed/"
#myDir2 = "/lorax/chenglab/yanding/Pub_Dat/GRch38/Fastq_files/"
myDir2 = "/lorax/chenglab/cc59/PubDat/organisms/human/genome/hg19/"

mytarDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/ATAC_seq_bin_sequence/"
dir.create(mytarDir1)

mychr = paste("chr", c(1:22, "X"),  sep="")
files = list.files(myDir1)
tag = grep(".Rda",files)
files = files[tag]
tag = grep("200bp",files)
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

[4]CNN model training

rm(list=ls())

#install.packages("tensorflow")
library(tensorflow)
#install_tensorflow()
library("Biostrings")
library("reticulate")

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/ATAC_seq_bin_sequence/broad_positive_200bp_ENCFF046YUJ_peak.fa"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/ATAC_seq_bin_sequence/broad_negative_200bp_ENCFF046YUJ_peak.fa"

library(abind)
library(keras)
#---------------------------
data <- readDNAStringSet(myinf1)
nn = length(data)/2
data = as.data.frame(data)
tmp1 = row.names(data)
tmp2 = data$x
names(tmp2) = tmp1
data = tmp2
res = array(0, c(length(data), 2000, 4))
for(k in 1:length(data))
{
	cat("\r", k)
	xx = unlist(strsplit(data[k], ""))
	se = which(xx=="A")
	res[k,se,1]=1
	se = which(xx=="T")
	res[k,se,2]=1
	se = which(xx=="C")
	res[k,se,3]=1
	se = which(xx=="G")
	res[k,se,4]=1	
}
dat.pos = res


data <- readDNAStringSet(myinf2)
tag = sample(1:length(data),20000)
data = data[tag]
nn = length(data)/2
data = as.data.frame(data)
tmp1 = row.names(data)
tmp2 = data$x
names(tmp2) = tmp1
data = tmp2
res = array(0, c(length(data), 2000, 4))
for(k in 1:length(data))
{
	cat("\r", k)
	xx = unlist(strsplit(data[k], ""))
	se = which(xx=="A")
	res[k,se,1]=1
	se = which(xx=="T")
	res[k,se,2]=1
	se = which(xx=="C")
	res[k,se,3]=1
	se = which(xx=="G")
	res[k,se,4]=1	
}
dat.neg = res



#---------------------------
nn1 = nrow(dat.pos)
se1 = sample(1:nn1)[1:(nn1/2)]
nn2 = nrow(dat.neg)
se2 = sample(1:nn2)[1:(nn2/2)]

tmp1 = dat.pos[se1, ,]
tmp2 = dat.neg[se2, ,]
tr.x = abind(tmp1, tmp2, along=1)
dim(tr.x)
tr.y = c(rep(1, dim(tmp1)[1]), rep(0, dim(tmp2)[1]))
idx = sample(1:length(tr.y))
tr.x = tr.x[idx,,]
tr.y = tr.y[idx]

tmp1 = dat.pos[-se1, ,]
tmp2 = dat.neg[-se2, ,]
te.x = abind(tmp1, tmp2, along=1)
dim(tr.x)
te.y = c(rep(1, dim(tmp1)[1]), rep(0, dim(tmp2)[1]))
idx = sample(1:length(tr.y))
te.x = te.x[idx,,]
te.y = te.y[idx]

#---------------------------
#  tr.y <- to_categorical(tr.y, 2)				## skip this for two-class classification
#tr.x <- array_reshape(tr.x, c(dim(tr.x)[1], 4, 100, 1))
#te.x <- array_reshape(te.x, c(dim(te.x)[1], 4, 100, 1))
dim(tr.x)
dim(te.x)

myoutf ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/Environment_for_CNN_200bp.Rda"
save.image(myoutf)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model <- keras_model_sequential() 
model %>% 
  layer_conv_1d(filter = 320, kernel_size = 8, activation = 'relu', input_shape = c(2000, 4)) %>% 
  layer_conv_1d(filter = 320, kernel_size = 8, activation = 'relu', input_shape = c(2000, 4)) %>% 
  layer_dropout(rate = 0.2) %>%
  layer_max_pooling_1d(pool_size = 4) %>%
  layer_conv_1d(filter = 480, kernel_size = 8, activation = 'relu', input_shape = c(2000, 4)) %>% 
  layer_conv_1d(filter = 480, kernel_size = 8, activation = 'relu', input_shape = c(2000, 4)) %>% 
  layer_dropout(rate = 0.2) %>%
  layer_max_pooling_1d(pool_size = 4) %>%
  layer_conv_1d(filter = 640, kernel_size = 8, activation = 'relu', input_shape = c(2000, 4)) %>%  
  layer_flatten() %>%
  layer_dropout(rate = 0.2) %>% 
  layer_dense(2000,activation = 'relu') %>% 
  layer_dense(120,activation = 'relu') %>% 
#  layer_dense(units = 2, activation = 'softmax')		## change units=1 and activation="sigmoid" for two-class classification
  layer_dense(units = 1, activation = 'sigmoid')

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#max pooling
#kernel size
#epoach
summary(model)

#set up checkpoints
checkpoint_path <- "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/CNN/cp_200bp.ckpt"

# Create checkpoint callback
cp_callback <- callback_model_checkpoint(
  filepath = checkpoint_path,
  monitor = "val_accuracy",
  save_weights_only = TRUE,
  save_best_only = TRUE,
  verbose = 1
)

model %>% compile(
#  loss = 'categorical_crossentropy',							## loss = "binary_crossentropy" for for two-class classification
  loss = 'binary_crossentropy',							
  optimizer = optimizer_rmsprop(lr = 0.0001, decay = 1e-6),
  metrics = c('accuracy')
)


history <- model %>% fit(
  tr.x, tr.y, 
  epochs = 50, batch_size = 128, 
  validation_split = 0.2,
  callbacks = list(cp_callback), # pass callback to training
)


myoutf ="/lorax/chenglab/yanding/ATAC_seq_integration/Figures/CNN/CNN_training.pdf"
pdf(width=5,height=5,myoutf)
plot(history)
dev.off()
###################################################################################################################
###################################################################################################################
###################################################################################################################

rm(list=ls())
library(liftOver)
library(tidyverse)
library(GenomicRanges)
library(ChIPpeakAnno)
library(rtracklayer)

###################################################
#define the convert function
peakDF2GRanges <- function(peak.df) {
    peak.gr=GRanges(seqnames=peak.df[,1],
        ranges=IRanges(peak.df[,2], peak.df[,3]))
    cn <- colnames(peak.df)
    if (length(cn) > 3) {
        for (i in 4:length(cn)) {
            mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
        }
    }
    return(peak.gr)
}

#define the liftover
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/Prostate/ATAC_seq/ENCFF046YUJ.bed"
data = read.table(myinf1,sep="\t",quote=NULL)
data = data[,c("V1","V2","V3")]
data = unique(data)
res1 = data
res1 = unique(res1)	
#liftover to hg19
res1 = peakDF2GRanges(res1)
seqlevelsStyle(res1) = "UCSC"
cur19 = liftOver(res1, ch)
cur19 = unlist(cur19)
genome(cur19) = "hg19"
cur19 = new("gwaswloc", cur19)
res1 = as.data.frame(cur19)
res1 = unique(res1)
	
myinf2 =  "/lorax/chenglab/yanding/Pub_Dat/Reference_genome/Hg19/Bin_200bp_hg19_genome.Rda"
load(myinf2)
tiles = as.data.frame(tiles)
tiles$start = tiles$start -1
tiles$end = tiles$end - 1

res1 = peakDF2GRanges(res1)
res2 = peakDF2GRanges(tiles)

library(ChIPpeakAnno)

library(GenomicRanges)

shared_peaks = findOverlapsOfPeaks(res1, res2, connectedPeaks="keepAll",minoverlap=100)
			
peak = shared_peaks$venn_cnt

all.peaks <- shared_peaks$all.peaks
res1.renamed <- all.peaks$res1
res2.renamed <- all.peaks$res2
res1 = as.data.frame(res1.renamed)
res2 = as.data.frame(res2.renamed)

peakNames <- melt(shared_peaks$peaklist[['res1///res2']]$peakNames, value.name="merged.peak.id")
peakNames <- peakNames$value.value

tag = which(row.names(res2) %in% peakNames)
res1 = res2[tag,]

tag = which(row.names(res2) %in% peakNames)
res2 = res2[-tag,]

res1 = as.data.frame(res1)
res2 = as.data.frame(res2)

res1 = res1[,c(1,2,3)]
res2 = res2[,c(1,2,3)]
colnames(res1) = colnames(res2) = c("V1","V2","V3")

#res1 = peakDF2GRanges(res1)
#res2 = peakDF2GRanges(res2)
#shared_peaks = findOverlapsOfPeaks(res1, res2, connectedPeaks="keepAll")
			
#peak = shared_peaks$venn_cnt

		
training_pos = as.data.frame(res1)
negative = as.data.frame(res2)


myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/Core_seq/core_positive_ENCFF046YUJ.Rda"
save(training_pos,file =myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/Core_seq/core_negative_ENCFF046YUJ.Rda"
save(negative,file =myoutf)

[2]Padding to 2000
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/Core_seq/core_positive_ENCFF046YUJ.Rda"
load(myinf1)

training_pos$start = as.numeric(as.vector(training_pos$V2)) - 900
training_pos$end = as.numeric(as.vector(training_pos$V3)) + 900
training_pos = training_pos[,c("V1","start","end")]
training_pos$ID = paste0(training_pos$V1,"_",training_pos$start,"_",training_pos$end)

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/Core_seq/core_negative_ENCFF046YUJ.Rda")

negative$start = as.numeric(as.vector(negative$V2)) - 900
negative$end = as.numeric(as.vector(negative$V3)) + 900
negative = negative[,c("V1","start","end" )]
negative$ID = paste0(negative$V1,"_",negative$start,"_",negative$end)

library(ChIPpeakAnno)

library(GenomicRanges)

peakDF2GRanges <- function(peak.df) {
    peak.gr=GRanges(seqnames=peak.df[,1],
        ranges=IRanges(peak.df[,2], peak.df[,3]))
    cn <- colnames(peak.df)
    if (length(cn) > 3) {
        for (i in 4:length(cn)) {
            mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
        }
    }
    return(peak.gr)
}


res1 = peakDF2GRanges(training_pos)
res2 = peakDF2GRanges(negative)
			
shared_peaks = findOverlapsOfPeaks(res1, res2, connectedPeaks="keepAll")
			
peak = shared_peaks$venn_cnt

peak
     res1 res2 Counts count.res1 count.res2
[1,]    0    0      0          0          0
[2,]    0    1  49607          0      49607
[3,]    1    0  19573      19573          0
[4,]    1    1    390        427        393

all.peaks <- shared_peaks$all.peaks
res1.renamed <- all.peaks$res1
res2.renamed <- all.peaks$res2
res1 = as.data.frame(res1.renamed)
res2 = as.data.frame(res2.renamed)

peakNames <- melt(shared_peaks$peaklist[['res1///res2']]$peakNames, value.name="merged.peak.id")
peakNames <- peakNames$value.value

tag = which(row.names(res2) %in% peakNames)
res2 = res2[-tag,]

res1 = as.data.frame(res1)
res2 = as.data.frame(res2)

res1 = res1[,c(1,2,3)]
res2 = res2[,c(1,2,3)]
colnames(res1) = colnames(res2) = c("V1","V2","V3")

tag = res1$V2 >0 & res1$V3 >0
res1 = res1[tag,]

tag = res2$V2 >0 & res2$V3 >0
res2 = res2[tag,]


training_pos = as.data.frame(res1)
negative = as.data.frame(res2)


myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/Padding/padding_positive_2000bp_ENCFF046YUJ.Rda"
save(training_pos,file =myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/Padding/padding_negative_2000bp_ENCFF046YUJ.Rda"
save(negative,file =myoutf)

[3]Get the sequence
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/Padding/"
#myDir2 = "/lorax/chenglab/yanding/Pub_Dat/GRch38/Fastq_files/"
myDir2 = "/lorax/chenglab/cc59/PubDat/organisms/human/genome/hg19/"

mytarDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/ATAC_seq_bin_sequence/"
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
	tag = sample(1:nrow(myp),20000)
	myp = myp[tag,]
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

[4]CNN Model

rm(list=ls())

#install.packages("tensorflow")
library(tensorflow)
#install_tensorflow()
library("Biostrings")
library("reticulate")

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/ATAC_seq_bin_sequence/padding_positive_2000bp_ENCFF046YUJ_peak.fa"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/ATAC_seq_bin_sequence/padding_negative_2000bp_ENCFF046YUJ_peak.fa"

library(abind)
library(keras)
#---------------------------
data <- readDNAStringSet(myinf1)
nn = length(data)/2
data = as.data.frame(data)
tmp1 = row.names(data)
tmp2 = data$x
names(tmp2) = tmp1
data = tmp2
res = array(0, c(length(data), 2000, 4))
for(k in 1:length(data))
{
	cat("\r", k)
	xx = unlist(strsplit(data[k], ""))
	se = which(xx=="A")
	res[k,se,1]=1
	se = which(xx=="T")
	res[k,se,2]=1
	se = which(xx=="C")
	res[k,se,3]=1
	se = which(xx=="G")
	res[k,se,4]=1	
}
dat.pos = res


data <- readDNAStringSet(myinf2)
nn = length(data)/2
data = as.data.frame(data)
tmp1 = row.names(data)
tmp2 = data$x
names(tmp2) = tmp1
data = tmp2
res = array(0, c(length(data), 2000, 4))
for(k in 1:length(data))
{
	cat("\r", k)
	xx = unlist(strsplit(data[k], ""))
	se = which(xx=="A")
	res[k,se,1]=1
	se = which(xx=="T")
	res[k,se,2]=1
	se = which(xx=="C")
	res[k,se,3]=1
	se = which(xx=="G")
	res[k,se,4]=1	
}
dat.neg = res



#---------------------------
nn1 = nrow(dat.pos)
se1 = sample(1:nn1)[1:(nn1/2)]
nn2 = nrow(dat.neg)
se2 = sample(1:nn2)[1:(nn2/2)]

tmp1 = dat.pos[se1, ,]
tmp2 = dat.neg[se2, ,]
tr.x = abind(tmp1, tmp2, along=1)
dim(tr.x)
tr.y = c(rep(1, dim(tmp1)[1]), rep(0, dim(tmp2)[1]))
idx = sample(1:length(tr.y))
tr.x = tr.x[idx,,]
tr.y = tr.y[idx]

tmp1 = dat.pos[-se1, ,]
tmp2 = dat.neg[-se2, ,]
te.x = abind(tmp1, tmp2, along=1)
dim(tr.x)
te.y = c(rep(1, dim(tmp1)[1]), rep(0, dim(tmp2)[1]))
idx = sample(1:length(tr.y))
te.x = te.x[idx,,]
te.y = te.y[idx]

#---------------------------
#  tr.y <- to_categorical(tr.y, 2)				## skip this for two-class classification
#tr.x <- array_reshape(tr.x, c(dim(tr.x)[1], 4, 100, 1))
#te.x <- array_reshape(te.x, c(dim(te.x)[1], 4, 100, 1))
dim(tr.x)
dim(te.x)

myoutf ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/Environment_for_CNN_200bp.Rda"
save.image(myoutf)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model <- keras_model_sequential() 
model %>% 
  layer_conv_1d(filter = 320, kernel_size = 10, activation = 'relu', input_shape = c(2000, 4)) %>% 
  layer_conv_1d(filter = 320, kernel_size = 10, activation = 'relu', input_shape = c(2000, 4)) %>% 
  layer_dropout(rate = 0.2) %>%
  layer_max_pooling_1d(pool_size = 2) %>%
  layer_conv_1d(filter = 480, kernel_size = 10, activation = 'relu', input_shape = c(2000, 4)) %>% 
  layer_conv_1d(filter = 480, kernel_size = 10, activation = 'relu', input_shape = c(2000, 4)) %>% 
  layer_dropout(rate = 0.2) %>%
  layer_max_pooling_1d(pool_size = 2) %>%
  layer_conv_1d(filter = 640, kernel_size = 8, activation = 'relu', input_shape = c(2000, 4)) %>%  
  layer_flatten() %>%
  layer_dropout(rate = 0.2) %>% 
  layer_dense(2000,activation = 'relu') %>% 
  layer_dense(120,activation = 'relu') %>% 
#  layer_dense(units = 2, activation = 'softmax')		## change units=1 and activation="sigmoid" for two-class classification
  layer_dense(units = 1, activation = 'sigmoid')



#max pooling
#kernel size
#epoach
summary(model)

#set up checkpoints
checkpoint_path <- "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Prostate/CNN/cp_200bp.ckpt"

# Create checkpoint callback
cp_callback <- callback_model_checkpoint(
  filepath = checkpoint_path,
  monitor = "val_accuracy",
  save_weights_only = TRUE,
  save_best_only = TRUE,
  verbose = 1
)

model %>% compile(
#  loss = 'categorical_crossentropy',							## loss = "binary_crossentropy" for for two-class classification
  loss = 'binary_crossentropy',							
  optimizer = optimizer_rmsprop(lr = 0.0001, decay = 1e-6),
  metrics = c('accuracy')
)


history <- model %>% fit(
  tr.x, tr.y, 
  epochs = 50, batch_size = 128, 
  validation_split = 0.2,
  callbacks = list(cp_callback), # pass callback to training
)


myoutf ="/lorax/chenglab/yanding/ATAC_seq_integration/Figures/CNN/CNN_training.pdf"
pdf(width=5,height=5,myoutf)
plot(history)
dev.off()