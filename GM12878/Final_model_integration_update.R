#####################
#Function define
#####################
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

[1]Histone 
library(randomForest)

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_neg.Rda"
load(myinf1)
pos_data = positive[validation_pos_tag,]
load(myinf2)
neg_data = negative[validation_neg_tag,]
data = rbind(pos_data, neg_data)
raw.data = data

data = raw.data
data = data[,2:ncol(data)]
pos = row.names(data)[seq(1,2000,1)]
neg = row.names(data)[seq(2001,4000,1)]
res = score_auc(data,pos,neg)

tmpinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/Histone_file_annotation_bigwig.txt"
info = read.table(tmpinf,sep="\t",quote=NULL,stringsAsFactors=F)

com = intersect(row.names(res),row.names(info))
info = info[com,]
res = res[com,]
info$AUC = res$AUC
tar = unique(info$Experiment.target)
raw.info = info

for(i in 1 : length(tar))
{
	cat("\r",i)
	info = raw.info
	tag = info$Experiment.target == tar[i]
	info = info[tag,]
	info = info[which.max(info$AUC),]
	
	if(i == 1)
	{
		final_info = info
	
	}else{
	
		final_info = rbind(final_info,info)
	}

}

data = raw.data
data = data[,c("tag",row.names(final_info))]

data$tag = factor(data$tag, levels = c("0","1"))
source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(data) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_integration_prediction.Rda"
save.image(file=myoutf)

#save the model for further usage
com = intersect(colnames(data),row.names(final_info))
data = data[,com]
final_info = final_info[com,]
colnames(data) = final_info$Experiment.target
tag = c(rep(1,2000),rep(0,2000))
data = cbind(tag,data)
data$tag = factor(data$tag,levels=c(0,1))

model <- randomForest(tag~., data=data, ntree=500,sampsize=c(2000,2000),importance=T)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_prediction_model.Rda"
save(model,file=myoutf)

[1] 0.7721785

[2]TF models
############################
#TF only
################################
library(randomForest)

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_TF_integration_Core.Rda")
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_neg.Rda"
load(myinf1)
pos_data = res2[validation_pos_tag,]
load(myinf2)
neg_data = res2[validation_neg_tag,]
data = rbind(pos_data, neg_data)
raw.data = data

data = raw.data
data = data[,2:ncol(data)]
pos = row.names(data)[seq(1,2000,1)]
neg = row.names(data)[seq(2001,4000,1)]
res = score_auc(data,pos,neg)

tmpinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/TF_file_annotation_bigWig.txt"
info = read.table(tmpinf,sep="\t",quote=NULL,stringsAsFactors=F)

com = intersect(row.names(res),row.names(info))
info = info[com,]
res = res[com,]
info$AUC = res$AUC
tar = unique(info$Experiment.target)
raw.info = info

for(i in 1 : length(tar))
{
	cat("\r",i)
	info = raw.info
	tag = info$Experiment.target == tar[i]
	info = info[tag,]
	info = info[which.max(info$AUC),]
	
	if(i == 1)
	{
		final_info = info
	
	}else{
	
		final_info = rbind(final_info,info)
	}

}

data = raw.data
data = data[,c("tag",row.names(final_info))]

data$tag = factor(data$tag, levels = c("0","1"))
source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(data) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda"
save.image(file=myoutf)

#save the model for further usage
com = intersect(colnames(data),row.names(final_info))
data = data[,com]
final_info = final_info[com,]
colnames(data) = final_info$Experiment.target
tag = c(rep(1,2000),rep(0,2000))
data = cbind(tag,data)
data$tag = factor(data$tag,levels=c(0,1))

model <- randomForest(tag~., data=data, ntree=500,sampsize=c(2000,2000),importance=T)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_prediction_model.Rda"
save(model,file=myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_prediction_data.Rda"
save(data,file=myoutf)


[3]Motif only
############################
#Motif only
################################
rm(list=ls())
mydir = "/lorax/chenglab/cc59/PubDat/organisms/human/motif/hsa_done/"
files = list.files(mydir)

library(randomForest)

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_neg.Rda")

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
colnames(res) = paste0("F_",seq(1,ncol(pos_count)))

data = cbind(tag,res)
data = as.data.frame(data)
data$tag = factor(data$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(data) 

#save the model for further usage
model <- randomForest(tag~., data=data, ntree=500,sampsize=c(2000,2000),importance=T)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Motif_integration_prediction.Rda"
save.image(file = myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Motif_prediction_model.Rda"
save(model,file=myoutf)

[4]Seq+TF models
rm(list=ls())
library(randomForest)
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/CNN/GM12878_local_CNN_model.Rda")
seq_data = prob
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda")
data = cbind(data,seq_data)

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(data) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Seq_TF_integration_prediction.Rda"
save.image(file=myoutf)

[5]Seq+HM models
rm(list=ls())
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/CNN/GM12878_local_CNN_model.Rda")
seq_data = prob

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_integration_prediction.Rda")
data = cbind(data,seq_data)

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(data) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Seq_Histone_integration_prediction.Rda"
save.image(file=myoutf)

[6]Seq+Motif models
rm(list=ls())
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/CNN/GM12878_local_CNN_model.Rda")
seq_data = prob

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Motif_integration_prediction.Rda")
data = cbind(data,seq_data)

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(data) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Motif_Seq_integration_prediction.Rda"
save.image(file=myoutf)

[7]Seq+Motif+HM+TF
rm(list=ls())
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/CNN/GM12878_local_CNN_model.Rda")
seq_data = prob

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Motif_integration_prediction.Rda")
motif_data = data[,2:ncol(data)]

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_integration_prediction.Rda")
histone_data = data[,2:ncol(data)]

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda")
tf_data = data[,2:ncol(data)]

data = cbind(seq_data,motif_data,histone_data,tf_data)
tag = c(rep(1,2000),rep(0,2000))
data = cbind(tag,data)
data$tag = factor(data$tag,levels=c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(data) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Seq_Motif_Histone_TF_integration_prediction.Rda"
save.image(file=myoutf)

[8]HM+TF
rm(list=ls())
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_integration_prediction.Rda")
histone_data = data[,2:ncol(data)]

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda")
tf_data = data[,2:ncol(data)]

data = cbind(histone_data,tf_data)
tag = c(rep(1,2000),rep(0,2000))
data = cbind(tag,data)
data$tag = factor(data$tag,levels=c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(data) 

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_TF_integration_prediction.Rda"
save.image(file=myoutf)




