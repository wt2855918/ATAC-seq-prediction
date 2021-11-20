[1]Create a binized genome
rm(list=ls())
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tiles <- tileGenome(seqinfo(txdb), tilewidth=100,
                    cut.last.tile.in.chrom=TRUE)
myoutf = "/lorax/chenglab/yanding/Pub_Dat/Reference_genome/Hg19/Bin_hg19_genome.Rda"
save(tiles, file = myoutf)

[2]Choose the random negative for GM12878
rm(list=ls())

f1 <- function(x, chr, sta, end) {paste(x[chr], x[sta]:x[end], sep="_")} #separate each base pair
f2 <- function(x, sig, sta, end) 
	{
		time = x[end]-x[sta]+1
		rep(as.character(x[sig]), time)
	}
bin.size = 100

myinf1 = "/lorax/chenglab/yanding/Pub_Dat/Reference_genome/Hg19/Bin_hg19_genome.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_positive.Rda"

load(myinf1)
tiles = as.data.frame(tiles)
tiles$start = tiles$start - 1
tiles$end = tiles$end -1
tiles$ID = paste0(tiles$seqnames,"_", tiles$start, "_",tiles$end )

load(myinf2)
positive$ID = paste0(positive$V1,"_", positive$V2,"_", positive$V3)

row.names(tiles) = tiles$ID
row.names(positive) = positive$ID

bin.size =100

#expand the positive by up and down stream 2kb
positive$V1 = as.vector(positive$V1)
positive$V2 = as.numeric(as.vector(positive$V2))
positive$V3 = as.numeric(as.vector(positive$V3))

positive$up_stream = positive$V2 - 2000
positive$dn_stream = positive$V3 + 2000

f2 <- function(x,sta,end){
	
	start = as.numeric(x[sta])
	end = as.numeric(x[end])
	
	seq(start,end,100)

}
com = intersect(row.names(tiles),row.names(positive))

tag = which(tiles$ID %in% com)
tiles = tiles[-tag,]

chr_name = unique(positive$V1)
tag = which(tiles$seqnames %in% chr_name)
tiles = tiles[tag,]

#tag = sample(seq(1,nrow(tiles)),50000)
#negative = tiles[tag,]

negative = tiles
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_random_negative.Rda"
save(negative, file = myoutf)

[3]Further generate the up and downstream for positive traning in GM12878
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda"
load(myinf1)
training_pos$up_stream = as.numeric(as.vector(training_pos$V2)) -2000
training_pos$dn_stream = as.numeric(as.vector(training_pos$V3)) + 2000

up_stream_start = training_pos$up_stream
up_stream_end = as.numeric(as.vector(training_pos$V2))
chrom = as.vector(training_pos$V1)
up_stream = mapply(function(x,y,z) paste0(x,"_",seq(y,z,100)), chrom, up_stream_start, up_stream_end)
up_stream = as.vector(up_stream)

dn_stream_end = training_pos$dn_stream
dn_stream_start = as.numeric(as.vector(training_pos$V3))
chrom = training_pos$V1
dn_stream = mapply(function(x,y,z) paste0(x,"_",seq(y,z,100)), chrom, dn_stream_start, dn_stream_end)
dn_stream = as.vector(dn_stream)

tar = c(up_stream, dn_stream)
tar = unique(tar)

tmpinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_random_negative.Rda"
load(tmpinf)
negative$ID_new = paste0(negative$seqnames, "_",negative$start)

tag = which(negative$ID_new %in% tar)
negative = negative[-tag,]

tag = sample(seq(1,nrow(negative)),50000)
negative = negative[tag,]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_random_negative_non_flank_overlap.Rda"
save(negative, file = myoutf)


[3]Generate the sequence for random negative GM12878
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/"
#myDir2 = "/lorax/chenglab/yanding/Pub_Dat/GRch38/Fastq_files/"
myDir2 = "/lorax/chenglab/cc59/PubDat/organisms/human/genome/hg19/"

mytarDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1/"
dir.create(mytarDir1)

mychr = paste("chr", c(1:22, "X","Y"),  sep="")
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
#k=3

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