rm(list=ls())
mydir = "/lorax/chenglab/cc59/PubDat/organisms/human/motif/hsa_done/"
files = list.files(mydir)

library(randomForest)
source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_neg.Rda")

bin.size = 100

pos = validation_pos
neg = validation_neg

#begin to calcualte the pos
pos$V2 = as.numeric(as.vector(pos$V2))
pos$V3 = as.numeric(as.vector(pos$V3))
	
pos$sta = ceiling((pos$V2+1)/bin.size)
pos$end = ceiling((pos$V3+1)/bin.size)
	
pos$ID = paste0(pos$V1,"_",pos$sta)
row.names(pos) = pos$ID
	
#begin to calulate the neg
	
neg$V2 = as.numeric(as.vector(neg$V2))
neg$V3 = as.numeric(as.vector(neg$V3))
	
neg$sta = ceiling((neg$V2+1)/bin.size)
neg$end = ceiling((neg$V3+1)/bin.size)
	
neg$ID = paste0(neg$V1,"_",neg$sta)
row.names(neg) = neg$ID
	
#create the matrix
motif_name = gsub(".txt","",files)

pos_count = matrix(0, nrow(pos), length(motif_name))	
row.names(pos_count) = row.names(pos)
colnames(pos_count) = motif_name
pos_count = as.data.frame(pos_count)

neg_count = matrix(0, nrow(neg), length(motif_name))
row.names(neg_count) = row.names(neg)
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

tag = c(rep(1,nrow(pos_count)),rep(0,nrow(neg_count)))
res = rbind(pos_count, neg_count)
colnames(res) = paste0("F_",seq(1,ncol(pos_count)))

data = cbind(tag,res)
data = as.data.frame(data)
data$tag = factor(data$tag, levels =c(0,1))
fit = ROC(data)

[1] 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_TF_motif_integration_AUC.Rda"
save.image(file = myoutf)
