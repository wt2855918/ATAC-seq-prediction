[1]Expand the positive region as training 
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda"
load(myinf1)
training_pos$start = as.numeric(as.vector(training_pos$V2)) - 950
training_pos$end = as.numeric(as.vector(training_pos$V3)) + 950
training_pos = training_pos[,c("V1","start","end")]

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_random_negative_non_flank_overlap.Rda")
negative$start = as.numeric(as.vector(negative$start)) - 950
negative$end = as.numeric(as.vector(negative$end)) + 950
negative = negative[,c("seqnames","start","end" )]

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

#Little leak which is fine 
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Broad_peak/broad_positive_ENCFF172DEA.Rda"
save(training_pos,file =myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Broad_peak/broad_negative_ENCFF172DEA.Rda"
save(negative,file =myoutf)

[2]Begin to get sequence
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Broad_peak/"
#myDir2 = "/lorax/chenglab/yanding/Pub_Dat/GRch38/Fastq_files/"
myDir2 = "/lorax/chenglab/cc59/PubDat/organisms/human/genome/hg19/"

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

[3]CNN training

rm(list=ls())

#install.packages("tensorflow")
library(tensorflow)
#install_tensorflow()
library("Biostrings")
library("reticulate")

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1/broad_positive_ENCFF172DEA_peak.fa"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1/broad_negative_ENCFF172DEA_peak.fa"

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

myoutf ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1/Environment_for_CNN.Rda"
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

model %>% compile(
#  loss = 'categorical_crossentropy',							## loss = "binary_crossentropy" for for two-class classification
  loss = 'binary_crossentropy',							
  optimizer = optimizer_rmsprop(lr = 0.0001, decay = 1e-6),
  metrics = c('accuracy')
)


history <- model %>% fit(
  tr.x, tr.y, 
  epochs = 100, batch_size = 128, 
  validation_split = 0.2
)

myoutf ="/lorax/chenglab/yanding/ATAC_seq_integration/Figures/CNN/CNN_training.pdf"
pdf(width=5,height=5,myoutf)
plot(history)
dev.off()

myoutf ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/CNN/GM12878_random_pading_1000bp_CNN.Rda"
save.image(myoutf)



