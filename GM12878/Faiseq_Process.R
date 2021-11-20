[5]Begin to create the control region
rm(list=ls())
#library(liftOver)
#library(tidyverse)
library(GenomicRanges)
library(ChIPpeakAnno)
library(rtracklayer)

myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/FAIRE-seq/"
myDir2 = "/lorax/chenglab/yanding/Pub_Dat/GRch38/Fastq_files/"

mytarDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/FAIRE_seq_peak_sequence_ctrl/"
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
	data = as.data.frame(data)
	data = unique(data)
	colnames(data) = c("V1","V2","V3","width")
	data = data[,c(1,2,3,4)]
	tag = data$width >100
	data = data[tag,]
	
	data$ID = paste0(data$V1,"_",data$V2,"_",data$V3)
	left_name = apply(data, 1, f1, chr="V1", sta="V2", end="V3", length = "width", upstream = 500)
	
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
	data = as.data.frame(data)
	data = unique(data)
	colnames(data) = c("V1","V2","V3","width")
	data = data[,c(1,2,3,4)]
	data = data[tag,]
	
	tag = data$width >100
	data = data[tag,]
	
	data$ID = paste0(data$V1,"_",data$V2,"_",data$V3)
	right_name = apply(data, 1, f1, chr="V1", sta="V2", end="V3", length = "width", dnstream = 500)
	
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

[2]Begin to examine the overlap
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


mydir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/FAIRE-seq/"
mydir2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/FAIRE_seq_peak_sequence_ctrl/"

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
	
	#liftover to hg19
	res1 = peakDF2GRanges(res1)
	res1 = as.data.frame(res1)
	res1 = unique(res1)
	
	tag = res1$width >100
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
	new_end = new_start + 99
	
	positive = cbind(chr_name, new_start, new_end)
	positive = as.data.frame(positive)
	colnames(positive) = c("V1","V2","V3")
	positive = unique(positive)
	
	myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/FAIRE_seq_peak_processed_bed/",my_nam[i],"_positive.Rda")
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
	new_end = new_start + 99
	
	negative = cbind(chr_name, new_start, new_end)
	negative = as.data.frame(negative)
	colnames(negative) = c("V1","V2","V3")
	negative = unique(negative)
	
	myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/FAIRE_seq_peak_processed_bed/",my_nam[i],"_negative.Rda")
	save(negative, file =myoutf)
	
}

[3]Begin to predict using Histone 
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_Histone_integration_Core.Rda"

load(myinf1)
res = res1
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/FAIRE_seq_peak_processed_bed/ENCFF001UYE_positive.Rda"

load(myinf2)
myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/FAIRE_seq_peak_processed_bed/ENCFF001UYE_negative.Rda"

load(myinf3)

positive$ID = paste0(positive$V1,"_",positive$V2,"_",positive$V3)
negative$ID = paste0(negative$V1,"_",negative$V2,"_",negative$V3)

sam1 = intersect(row.names(res),positive$ID)
sam2 = intersect(row.names(res),negative$ID)

res1 = res[sam1,]
res2 = res[sam2,]

#tag = res1$tag == "1"
#res1 = res1[tag==0,]

res = rbind(res1,res2)

raw.res1 = res1
raw.res2 = res2

#examine the share histone
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_integration_prediction.Rda")

tag1 = sample(seq(1,nrow(raw.res1)),2000)
tag2 = sample(seq(1,nrow(raw.res2)),2000)

res = rbind(raw.res1[tag1,],raw.res2[tag2,])

com = intersect(colnames(data),colnames(res))
res = res[,com]

tag = c(rep(1,2000),rep(0,2000))
res$tag = tag
res$tag = factor(res$tag, levels=c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

#using the self prediction
HM_fit = ROC(res) 

#AUC curve cell line
label <- paste0("AUC=",round(mean(1-fit[[3]]),2))
color <- brewer.pal(10,"Set1")[1]
myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/Histone_FAIR_prediction.pdf"
pdf(myoutf, width= 5, height= 5)
par(pty="s")
plot(1-fit[[1]],1-fit[[2]],main="Histone", xlab="False Positive Rate (1-Specificity)",ylab="True Positive Rate (Sensitivity)",
xlim=c(0,1),ylim=c(0,1),cex=0,cex.main=1.5,cex.lab=1,font.main=2) #xaxs="i",yaxs="i"
par(new=TRUE)
lines(1-fit[[1]],1-fit[[2]],lwd=2,col=color[1],lty=1)
abline(0,1,lty=3)
legend("bottomright",label,lty=rep(1,1),lwd=rep(2.5,1),col=color[1],cex=0.75, box.lty=0)
dev.off()

#change the name 
com = intersect(colnames(res), row.names(final_info))
res = res[,com]
final_info = final_info[com,]
colnames(res) = final_info$Experiment.target
tag = c(rep(1,2000),rep(0,2000))
res = cbind(tag,res)
res$tag = factor(res$tag, levels=c(0,1))


#using the model prediction
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_prediction_model.Rda")
data = res
res = NULL
tmp = predict(model, data[,-1], type="prob")
tmp = data.frame(data[,1], tmp)
res = rbind(res, tmp)

thr = (1:99)*0.01
yy =  xx =  rep(0, length(thr))
fdr = rep(0,99)
for(i in 1:length(thr))
{
	aa = sum(res[,2]>=thr[i] & res[,1]=="1")
	bb = sum(res[,2]<thr[i] & res[,1]=="1" )
	cc = sum(res[,2]>=thr[i] & res[,1]=="0")
	dd = sum(res[,2]<thr[i] & res[,1]=="0")
	fdr[i] = aa/sum(res[,2]>=thr[i])
	yy[i] = aa/(aa+bb)
	xx[i] = cc/(cc+dd)
}
xx = c(1, xx, 0) #TPR
yy = c(1, yy, 0) #FPR
tmp1 = tmp2 = rep(0,100)
for(i in 1:100)
{
	tmp1[i] = xx[i]-xx[i+1]
	tmp2[i] = (yy[i+1]+yy[i])/2	
}
myauc = sum(tmp1*tmp2)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_FAIR_prediction_AUC.Rda"
save.image(myoutf)


[4]Begin to predict using TF
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_TF_integration_Core.Rda"

load(myinf1)
res = res2

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/FAIRE_seq_peak_processed_bed/ENCFF001UYE_positive.Rda"
load(myinf2)

myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/FAIRE_seq_peak_processed_bed/ENCFF001UYE_negative.Rda"
load(myinf3)

positive$ID = paste0(positive$V1,"_",positive$V2,"_",positive$V3)
negative$ID = paste0(negative$V1,"_",negative$V2,"_",negative$V3)

sam1 = intersect(row.names(res),positive$ID)
sam2 = intersect(row.names(res),negative$ID)

res1 = res[sam1,]
res2 = res[sam2,]

#tag = res1$tag == "1"
#res1 = res1[tag==0,]

res = rbind(res1,res2)

raw.res1 = res1
raw.res2 = res2

#examine the share histone
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda")

tag1 = sample(seq(1,nrow(raw.res1)),2000)
tag2 = sample(seq(1,nrow(raw.res2)),2000)

res = rbind(raw.res1[tag1,],raw.res2[tag2,])

com = intersect(colnames(data),colnames(res))
res = res[,com]

tag = c(rep(1,2000),rep(0,2000))
res$tag = tag
res$tag = factor(res$tag, levels=c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

TF_fit = ROC(res) 

#AUC curve cell line
label <- paste0("AUC=",round(mean(1-fit[[3]]),2))
color <- brewer.pal(10,"Set1")[1]
myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/TF_FAIR_prediction.pdf"
pdf(myoutf, width= 5, height= 5)
par(pty="s")
plot(1-fit[[1]],1-fit[[2]],main="Histone", xlab="False Positive Rate (1-Specificity)",ylab="True Positive Rate (Sensitivity)",
xlim=c(0,1),ylim=c(0,1),cex=0,cex.main=1.5,cex.lab=1,font.main=2) #xaxs="i",yaxs="i"
par(new=TRUE)
lines(1-fit[[1]],1-fit[[2]],lwd=2,col=color[1],lty=1)
abline(0,1,lty=3)
legend("bottomright",label,lty=rep(1,1),lwd=rep(2.5,1),col=color[1],cex=0.75, box.lty=0)
dev.off()

#change the name 
com = intersect(colnames(res), row.names(final_info))
res = res[,com]
final_info = final_info[com,]
colnames(res) = final_info$Experiment.target
tag = c(rep(1,2000),rep(0,2000))
res = cbind(tag,res)
res$tag = factor(res$tag, levels=c(0,1))

#using the model prediction
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_prediction_model.Rda")
data = res
res = NULL
tmp = predict(model, data[,-1], type="prob")
tmp = data.frame(data[,1], tmp)
res = rbind(res, tmp)

thr = (1:99)*0.01
yy =  xx =  rep(0, length(thr))
fdr = rep(0,99)
for(i in 1:length(thr))
{
	aa = sum(res[,2]>=thr[i] & res[,1]=="1")
	bb = sum(res[,2]<thr[i] & res[,1]=="1" )
	cc = sum(res[,2]>=thr[i] & res[,1]=="0")
	dd = sum(res[,2]<thr[i] & res[,1]=="0")
	fdr[i] = aa/sum(res[,2]>=thr[i])
	yy[i] = aa/(aa+bb)
	xx[i] = cc/(cc+dd)
}
xx = c(1, xx, 0) #TPR
yy = c(1, yy, 0) #FPR
tmp1 = tmp2 = rep(0,100)
for(i in 1:100)
{
	tmp1[i] = xx[i]-xx[i+1]
	tmp2[i] = (yy[i+1]+yy[i])/2	
}
myauc = sum(tmp1*tmp2)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_FAIR_prediction_AUC.Rda"
save.image(myoutf)

