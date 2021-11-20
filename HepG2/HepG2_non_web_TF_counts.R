[1]Examine the positive not included in the ENCODE annotation
rm(list=ls())
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
files = list.files(myinf)
tag = grep(".tab",files)
files = files[tag]
file_nam = gsub(".tab","",files)

tmpinf =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/TF_file_annotation_bigWig.txt"
info = read.table(tmpinf,sep="\t",quote=NULL)

tag = which(file_nam %in% row.names(info))
files = files[-tag]
file_nam = file_nam[-tag]

#load the positive and negative
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/HepG2_training_and_validation_pos.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/HepG2_training_and_validation_neg.Rda"
	
load(myinf1)
load(myinf2)
	
all_bins = c(validation_pos_tag, validation_neg_tag)
res = matrix(0,length(all_bins),length(files))
row.names(res) = all_bins
colnames(res) = file_nam
res = as.data.frame(res)

myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"

myk = kkkk

tmpinf = paste0(myinf,files[myk])
data = read.table(tmpinf,sep="\t",quote=NULL)
data[,"V3"] = data[,"V3"] -1
	
row.names(data) = paste0(data[,"V1"],"_",data[,"V2"],"_",data[,"V3"])	
data = data[all_bins,]
	
myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_counts_bins_logFC_non_web/",file_nam[myk],".Rda")
save(data,file = myoutf)

[2]Begin to modifiy

rm(list=ls())
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
files = list.files(myinf)
tag = grep(".tab",files)
files = files[tag]
file_nam = gsub(".tab","",files)

tmpinf =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/TF_file_annotation_bigWig.txt"
info = read.table(tmpinf,sep="\t",quote=NULL)

tag = which(file_nam %in% row.names(info))
files = files[-tag]
file_nam = file_nam[-tag]

fnum = length(file_nam)

myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/HepG2_TF_non_web/"
dir.create(myDir1)
mytemplate = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Template/HepG2/Tempalte_for_counting_non_web_HepG2_TF_sig.R"
mysub = paste(myDir1, "submit.sp", sep="")
setwd(myDir1)

conIn = file(mytemplate, "r")
rawdata = readLines(conIn)
close(conIn)

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

[3]Begin to merge
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_counts_bins_logFC_non_web/"
files = list.files(mydir)

my_nam = gsub(".Rda","",files)
res = matrix(0, 4000, length(files))
colnames(res) = my_nam
res = as.data.frame(res)

for(i in 1 : length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir,files[i])
	load(tmpinf)
	
    res[,i] = data$V4
}

#define function

score_auc <-function(data,pos, neg){

	res = matrix(0, ncol(data), 8)
	colnames(res) = c("pCR.Avg", "RD.Avg", "Tscore", "pval.t", "pval.w", "FDR.w", "AUC.ori", "AUC")
	row.names(res) = colnames(data)
	res = as.data.frame(res)
	dat1 = data[row.names(data)%in%pos,]
	dat2 = data[row.names(data)%in%neg,]
	res[,1] = apply(dat1, 2, function(x) mean(x,na.rm=T))
	res[,2] = apply(dat2, 2, function(x) mean(x,na.rm=T))

	for(k in 1:ncol(data))
	{
		cat("\r",k)
		tmp = t.test(dat1[, k], dat2[,k])
		res[k,3] = tmp$statistic
		res[k,4] = tmp$p.value
		tmp = wilcox.test(dat1[, k], dat2[,k])
		res[k,5] = tmp$p.value

		xx = data[,k]
		names(xx) = row.names(data)
		xx = sort(xx)
		xx= names(xx)%in%pos
		fp = 1-xx 
		tp = xx
		for(j in length(xx):2)
		{
			fp[j-1]= fp[j]+fp[j-1]
			tp[j-1]= tp[j]+tp[j-1]
		}
		fp = fp/length(neg)
		tp = tp/length(pos)	
		xx = c(1, fp, 0)
		yy = c(1, tp, 0)
		tmp1 = tmp2 = rep(0,length(xx)-1)
		for(i in 1:length(tmp1))
		{
			tmp1[i] = xx[i]-xx[i+1]
			tmp2[i] = (yy[i+1]+yy[i])/2	
		}
		res[k,7] = sum(tmp1*tmp2)
	}
	res[,6] = p.adjust(res[,5], method="BH")
	res[,8] = ifelse(res[,7]>0.5, res[,7], 1-res[,7])

	return(res)
}

pos = row.names(res)[1:2000]
neg = row.names(res)[2001:4000]

info = score_auc(res,pos,neg)

[4]Examine the file accession difference
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/metadata.tsv"
res1 = read.table(myinf1,sep="\t",quote=NULL,header=T)

library(ENCODExplorer)
res <- queryEncode(organism = "Homo sapiens", assay = "ChIP-seq",
                      biosample_name = "HepG2", file_format = "bigWig",
                      fixed = TRUE)                     
tag = res$output_type  == "fold change over control"
res = res[tag,]
tag = grep("hg19",res$assembly)
res = res[tag,]
tag = grep("transcription factor",res$investigated_as)
res = res[tag,]
res2 = res

row.names(res1) = res1$file_accession
row.names(res2) = res2$file_accession








