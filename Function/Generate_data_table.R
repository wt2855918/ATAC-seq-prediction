
[1]Motif GM12878
rm(list=ls())
mydir = "/mount/ictr1/chenglab/cc59/PubDat/organisms/human/motif/hsa_done/"
files = list.files(mydir)

library(randomForest)

load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda")
load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_neg.Rda")

bin.size = 100

pos = validation_pos_tag
neg = validation_neg_tag

pos = as.vector(pos)
neg = as.vector(neg)

#begin to calcualte the pos
chr = lapply(pos, function(x) strsplit(x,"_")[[1]][1])
start = lapply(pos, function(x) strsplit(x,"_")[[1]][2])

chr = as.vector(unlist(chr))
start = as.numeric(as.vector(unlist(start)))
new_start = ceiling((start+1)/bin.size)
pos_ID = paste0(chr,"_",new_start)
pos_nam = pos

#begin to calulate the neg
chr = lapply(neg, function(x) strsplit(x,"_")[[1]][1])
start = lapply(neg, function(x) strsplit(x,"_")[[1]][2])

chr = as.vector(unlist(chr))
start = as.numeric(as.vector(unlist(start)))
new_start = ceiling((start+1)/bin.size)
neg_ID = paste0(chr,"_",new_start)
neg_nam = neg

#create the matrix
motif_name = gsub(".txt","",files)

pos_count = matrix(0, length(pos_ID), length(motif_name))	
row.names(pos_count) = pos_ID
colnames(pos_count) = motif_name
pos_count = as.data.frame(pos_count)

neg_count = matrix(0, length(neg_ID), length(motif_name))
row.names(neg_count) = neg_ID
colnames(neg_count) = motif_name
neg_count = as.data.frame(neg_count)
	
for(i in 1 : length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir,files[i])
	tmp = read.table(tmpinf,sep="\t",quote=NULL)
	
	pos_com = intersect(row.names(pos_count),tmp$V1)
	neg_com = intersect(row.names(neg_count),tmp$V1)
	
	pos_count[pos_com,i]= 1
	neg_count[neg_com,i]= 1
	
}

#recover the names
row.names(pos_count) = pos_nam
row.names(neg_count) = neg_nam


tag = c(rep(1,nrow(pos_count)),rep(0,nrow(neg_count)))
res = rbind(pos_count, neg_count)
#colnames(res) = paste0("F_",seq(1,ncol(pos_count)))

data = cbind(tag,res)
data = as.data.frame(data)
data$tag = factor(data$tag, levels =c(0,1))

myoutf = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/Finished_Table/GM12878/GM12878_motif_integration.Rda"
save(data,file = myoutf)

[2]Motif HepG2
rm(list=ls())
mydir = "/mount/ictr1/chenglab/cc59/PubDat/organisms/human/motif/hsa_done/"
files = list.files(mydir)

library(randomForest)

load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/HepG2/HepG2_training_and_validation_pos.Rda")
load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/HepG2/HepG2_training_and_validation_neg.Rda")

bin.size = 100

pos = validation_pos_tag
neg = validation_neg_tag

pos = as.vector(pos)
neg = as.vector(neg)

#begin to calcualte the pos
chr = lapply(pos, function(x) strsplit(x,"_")[[1]][1])
start = lapply(pos, function(x) strsplit(x,"_")[[1]][2])

chr = as.vector(unlist(chr))
start = as.numeric(as.vector(unlist(start)))
new_start = ceiling((start+1)/bin.size)
pos_ID = paste0(chr,"_",new_start)
pos_nam = pos

#begin to calulate the neg
chr = lapply(neg, function(x) strsplit(x,"_")[[1]][1])
start = lapply(neg, function(x) strsplit(x,"_")[[1]][2])

chr = as.vector(unlist(chr))
start = as.numeric(as.vector(unlist(start)))
new_start = ceiling((start+1)/bin.size)
neg_ID = paste0(chr,"_",new_start)
neg_nam = neg

#create the matrix
motif_name = gsub(".txt","",files)

pos_count = matrix(0, length(pos_ID), length(motif_name))	
row.names(pos_count) = pos_ID
colnames(pos_count) = motif_name
pos_count = as.data.frame(pos_count)

neg_count = matrix(0, length(neg_ID), length(motif_name))
row.names(neg_count) = neg_ID
colnames(neg_count) = motif_name
neg_count = as.data.frame(neg_count)
	
for(i in 1 : length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir,files[i])
	tmp = read.table(tmpinf,sep="\t",quote=NULL)
	
	pos_com = intersect(row.names(pos_count),tmp$V1)
	neg_com = intersect(row.names(neg_count),tmp$V1)
	
	pos_count[pos_com,i]= 1
	neg_count[neg_com,i]= 1
	
}

#recover the names
row.names(pos_count) = pos_nam
row.names(neg_count) = neg_nam


tag = c(rep(1,nrow(pos_count)),rep(0,nrow(neg_count)))
res = rbind(pos_count, neg_count)
#colnames(res) = paste0("F_",seq(1,ncol(pos_count)))

data = cbind(tag,res)
data = as.data.frame(data)
data$tag = factor(data$tag, levels =c(0,1))

myoutf = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/Finished_Table/HepG2/HepG2_motif_integration.Rda"
save(data,file = myoutf)

[3]GM12878 TF and histone and sequence
rm(list=ls())
load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/HepG2/CNN/HepG2_local_CNN_model.Rda")
seq_data = prob

load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/Finished_Table/HepG2/HepG2_motif_integration.Rda")
motif_data = data[,2:ncol(data)]

load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_integration_prediction.Rda")
histone_data = data[,2:ncol(data)]

tmpinf1 = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/Annotation/HepG2/Histone_file_annotation_bigwig.txt"
res1 = read.table(tmpinf1,sep="\t",quote=NULL)
#target = unique(res1$Experiment.target)
#for(i in 1 : length(target))
#{
#	cat("\r",i)
#	tag = res1$Experiment.target == target[i]
#	name = paste0(target[i],"_Rep_",seq(1,sum(tag),1))
#	res1$Experiment.target[tag] = name
#}

row.names(res1) = res1$File.accession
com = intersect(colnames(histone_data),row.names(res1))
histone_data = histone_data[,com]
res1 = res1[com,]
colnames(histone_data) = res1$Experiment.target

load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/HepG2/Final_ID/TF_integration_prediction.Rda")
tf_data = data[,2:ncol(data)]

tmpinf2 =  "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/Annotation/HepG2/TF_file_annotation_bigWig.txt"
res2 = read.table(tmpinf2,sep="\t",quote=NULL)
#target = unique(res2$Experiment.target)
#for(i in 1 : length(target))
#{
#	cat("\r",i)
#	tag = res2$Experiment.target == target[i]
#	name = paste0(target[i],"_Rep_",seq(1,sum(tag),1))
#	res2$Experiment.target[tag] = name
#}

row.names(res2) = res2$File.accession
com = intersect(colnames(tf_data),row.names(res2))
tf_data = tf_data[,com]
res2 = res2[com,]
colnames(tf_data) = res2$Experiment.target

data = cbind(seq_data,motif_data,histone_data,tf_data)

myoutf = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/HepG2/Final_ID/Final_Integration_table_HepG2.txt"
write.table(data,myoutf,sep="\t",quote=F)


[2]GM12878
rm(list=ls())
load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/GM12878/CNN/GM12878_local_CNN_model.Rda")
seq_data = prob

load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/Finished_Table/GM12878/GM12878_motif_integration.Rda")
motif_data = data[,2:ncol(data)]

load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_integration_prediction.Rda")
histone_data = data[,2:ncol(data)]

tmpinf1 = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/Annotation/GM12878/Histone_file_annotation_bigwig.txt"
res1 = read.table(tmpinf1,sep="\t",quote=NULL)
#target = unique(res1$Experiment.target)
#for(i in 1 : length(target))
#{
#	cat("\r",i)
#	tag = res1$Experiment.target == target[i]
#	name = paste0(target[i],"_Rep_",seq(1,sum(tag),1))
#	res1$Experiment.target[tag] = name
#}

row.names(res1) = res1$File.accession
com = intersect(colnames(histone_data),row.names(res1))
histone_data = histone_data[,com]
res1 = res1[com,]
colnames(histone_data) = res1$Experiment.target

load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda")
tf_data = data[,2:ncol(data)]

tmpinf2 =  "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/Annotation/GM12878/TF_file_annotation_bigWig.txt"
res2 = read.table(tmpinf2,sep="\t",quote=NULL)
#target = unique(res2$Experiment.target)
#for(i in 1 : length(target))
#{
#	cat("\r",i)
#	tag = res2$Experiment.target == target[i]
#	name = paste0(target[i],"_Rep_",seq(1,sum(tag),1))
#	res2$Experiment.target[tag] = name
#}

row.names(res2) = res2$File.accession
com = intersect(colnames(tf_data),row.names(res2))
tf_data = tf_data[,com]
res2 = res2[com,]
colnames(tf_data) = res2$Experiment.target

data = cbind(seq_data,motif_data,histone_data,tf_data)

myoutf = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/HepG2/Final_ID/Final_Integration_table_GM12878.txt"
write.table(data,myoutf,sep="\t",quote=F)

[3]Begin to calculate enrichment analysis GM12878
rm(list=ls())
load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_integration_prediction.Rda")
#histone_data = data[,2:ncol(data)]

res = matrix(0,(ncol(data)-1),4)
colnames(res) = c("Avg_open","Avg_close","T.score","P.value")
row.names(res) = colnames(data)[2:ncol(data)]

for(i in 2 : ncol(data))
{
	tag = data$tag == 1
	
	xx1 = data[tag==1,i]
	xx2 = data[tag==0,i]

	fit = t.test(xx1,xx2)
	
	res[i-1,"Avg_open"] = mean(xx1)
	res[i-1,"Avg_close"] = mean(xx2)
	res[i-1,"T.score"] = fit$statistic
	res[i-1,"P.value"] = fit$p.value
}

[3]Begin to calculate enrichment analysis GM12878
rm(list=ls())
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

library(randomForest)

load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/HepG2/Integration/HepG2_TF_integration_Core.Rda")
myinf1 = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/HepG2/HepG2_training_and_validation_pos.Rda"
myinf2 = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/HepG2/HepG2_training_and_validation_neg.Rda"
load(myinf1)
pos_data = res2[validation_pos_tag,]
load(myinf2)
neg_data = res2[validation_neg_tag,]
data = rbind(pos_data, neg_data)
raw.data = data

data = raw.data
data = data[,2:ncol(data)]
data = apply(data,2,function(x) (x-mean(x))/sd(x))

pos = row.names(data)[seq(1,2000,1)]
neg = row.names(data)[seq(2001,4000,1)]
res = score_auc(data,pos,neg)
res = res[,1:2]

xx1 = quantile(as.matrix(res),0.9)
xx2 = quantile(as.matrix(res),0.1)
tag1 = res >= xx1
tag2 = res <= xx2
res[tag1] = xx1
res[tag2] = xx2

myoutf = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Figures/Suppl_Figure/HepG2_TF_heatmap.pdf"
pdf(myoutf,width=3,height=5)
pheatmap(res,show_rownames = F, show_colnames =F)
dev.off()

[2]HepG2 Histone
library(randomForest)

myinf1 = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/HepG2/HepG2_training_and_validation_pos.Rda"
myinf2 = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/HepG2/HepG2_training_and_validation_neg.Rda"
load(myinf1)
pos_data = positive[validation_pos_tag,]
load(myinf2)
neg_data = negative[validation_neg_tag,]
data = rbind(pos_data, neg_data)
raw.data = data

data = raw.data
data = data[,2:ncol(data)]
data = apply(data,2,function(x) (x-mean(x))/sd(x))

pos = row.names(data)[seq(1,2000,1)]
neg = row.names(data)[seq(2001,4000,1)]
res = score_auc(data,pos,neg)
res = res[,1:2]

xx1 = quantile(as.matrix(res),0.9)
xx2 = quantile(as.matrix(res),0.1)
tag1 = res >= xx1
tag2 = res <= xx2
res[tag1] = xx1
res[tag2] = xx2

myoutf = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Figures/Suppl_Figure/HepG2_Histone_TF_heatmap.pdf"
pdf(myoutf,width=3,height=5)
pheatmap(res,show_rownames = F, show_colnames =F)
dev.off()

[3]HepG2 motif enrichment
tmpinf = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/Finished_Table/HepG2/HepG2_motif_integration.Rda"
load(tmpinf)

data = data[,2:ncol(data)]
res1 = data[1:2000,]
res2 = data[2001:4000,]

xx1 = apply(res1,2,sum)
xx2 = apply(res2,2,sum)

ER = (xx1/2000)/((xx1+xx2)/4000)
res = as.data.frame(ER)
res$ER_1 = (xx2/2000)/((xx1+xx2)/4000)
res = res[order(res$ER,decreasing=T),]

myoutf = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Figures/Suppl_Figure/HepG2_Motif_heatmap.pdf"
pdf(myoutf,width=3,height=5)
pheatmap(res,show_rownames = F, show_colnames =F, cluster_rows=F)
dev.off()

[4]GM12878 TF
library(randomForest)

load("/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/GM12878/Integration/GM12878_TF_integration_Core.Rda")
myinf1 = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda"
myinf2 = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_neg.Rda"
load(myinf1)
pos_data = res2[validation_pos_tag,]
load(myinf2)
neg_data = res2[validation_neg_tag,]
data = rbind(pos_data, neg_data)
raw.data = data

data = raw.data
data = data[,2:ncol(data)]
data = apply(data,2,function(x) (x-mean(x))/sd(x))

avg = apply(data,2,mean)
tag = is.na(avg)
data = data[,tag==0]

pos = row.names(data)[seq(1,2000,1)]
neg = row.names(data)[seq(2001,4000,1)]
res = score_auc(data,pos,neg)
res = res[,1:2]

xx1 = quantile(as.matrix(res),0.9)
xx2 = quantile(as.matrix(res),0.1)
tag1 = res >= xx1
tag2 = res <= xx2
res[tag1] = xx1
res[tag2] = xx2

myoutf = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Figures/Suppl_Figure/GM12878_TF_heatmap.pdf"
pdf(myoutf,width=3,height=5)
pheatmap(res,show_rownames = F, show_colnames =F)
dev.off()

[5]GM12878 Histone 
library(randomForest)

myinf1 = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda"
myinf2 = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_neg.Rda"
load(myinf1)
pos_data = positive[validation_pos_tag,]
load(myinf2)
neg_data = negative[validation_neg_tag,]
data = rbind(pos_data, neg_data)
raw.data = data

data = raw.data
data = data[,2:ncol(data)]
data = apply(data,2,function(x) (x-mean(x))/sd(x))

pos = row.names(data)[seq(1,2000,1)]
neg = row.names(data)[seq(2001,4000,1)]
res = score_auc(data,pos,neg)
res = res[,1:2]

xx1 = quantile(as.matrix(res),0.9)
xx2 = quantile(as.matrix(res),0.1)
tag1 = res >= xx1
tag2 = res <= xx2
res[tag1] = xx1
res[tag2] = xx2

myoutf = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Figures/Suppl_Figure/GM12878_Histone_TF_heatmap.pdf"
pdf(myoutf,width=3,height=5)
pheatmap(res,show_rownames = F, show_colnames =F)
dev.off()

[6]GM12878 Motif 
tmpinf = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Data/Finished_Table/GM12878/GM12878_motif_integration.Rda"
load(tmpinf)

data = data[,2:ncol(data)]
res1 = data[1:2000,]
res2 = data[2001:4000,]

xx1 = apply(res1,2,sum)
xx2 = apply(res2,2,sum)

ER = (xx1/2000)/((xx1+xx2)/4000)
res = as.data.frame(ER)
res$ER_1 = rep(0,nrow(res))
res = res[order(res$ER,decreasing=T),]

myoutf = "/mount/cheng-scratch/yanding/Dart/ATAC_seq_integration/Figures/Suppl_Figure/GM12878_Motif_heatmap.pdf"
pdf(myoutf,width=3,height=5)
pheatmap(res,show_rownames = F, show_colnames =F, cluster_rows=F)
dev.off()

