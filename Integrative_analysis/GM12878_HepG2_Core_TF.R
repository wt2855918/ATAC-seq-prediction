[1]GM12878
rm(list=ls())

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_KO_RF/"
files = list.files(mydir)

tmpinf =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_TF_fit.Rda"

my_nam = gsub("_KO_TF.Rda","",files)

sample_type = sapply(my_nam, function(x) strsplit(x,"_")[[1]][1])
Order_type = sapply(my_nam, function(x) strsplit(x,"_")[[1]][3])
sample_name = sapply(my_nam, function(x) strsplit(x,"_")[[1]][4])

info = as.data.frame(cbind(sample_type, Order_type, sample_name))
info$AUC = rep(0,nrow(info))

for(i in 1 : length(files))
{
	cat("\r",i)
	myinf = paste0(mydir, files[i])
	load(myinf)
	
	if(mean(fit[[3]])< 0.5)
	{
		AUC = 1- mean(fit[[3]])

	}else{
	
		AUC = mean(fit[[3]])
	}
	
	info$AUC[i] = AUC
}

info$Target = info$sample_name 
KO_info = info

load(tmpinf)
all_target = c("GM12878","0","All",mean(1-fit[[3]]),"All")

info = KO_info
info = apply(info,2, function(x) as.vector(x))
info = rbind(all_target, info)
info = as.data.frame(info)

row.names(info) = info$sample_name
info$Order_type = as.numeric(as.vector(info$Order_type))
info$AUC = as.numeric(as.vector(info$AUC))

info = info[order(info$Order_type),]
#info$Name = paste0(info$sample_name,"_",info$Target)
info$Name = info$Target
info$Name = factor(info$Name, levels = info$Name)
info$AUC = round(info$AUC,2)
AUC_info = info

library(randomForest)
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_TF_data.Rda"
load(myinf1)

label = GM12878_data$tag


tag = AUC_info$AUC < 0.80
target = AUC_info$Target[tag==0]
target = as.vector(target)
tag = grep("All",target)
target = target[-tag]

tag = which(colnames(GM12878_data) %in% target)
GM12878_data = GM12878_data[,-tag]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_core_TF_data.Rda"
save(GM12878_data, file = myoutf)

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(GM12878_data) 
GM12878_fit = fit

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_core_TF_fit.Rda"
save(GM12878_fit,file = myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_core_TF_model.Rda"
GM12878_model <- randomForest(tag~., data=GM12878_data, ntree=500,sampsize=c(2000,2000),importance=T)
save(GM12878_model,file = myoutf)

[2]HepG2
rm(list=ls())

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_KO_RF/"
files = list.files(mydir)

tmpinf =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/HepG2_shared_TF_fit.Rda"

my_nam = gsub("_KO_TF.Rda","",files)

sample_type = sapply(my_nam, function(x) strsplit(x,"_")[[1]][1])
Order_type = sapply(my_nam, function(x) strsplit(x,"_")[[1]][3])
sample_name = sapply(my_nam, function(x) strsplit(x,"_")[[1]][4])

info = as.data.frame(cbind(sample_type, Order_type, sample_name))
info$AUC = rep(0,nrow(info))

for(i in 1 : length(files))
{
	cat("\r",i)
	myinf = paste0(mydir, files[i])
	load(myinf)
	
	if(mean(fit[[3]])< 0.5)
	{
		AUC = 1- mean(fit[[3]])

	}else{
	
		AUC = mean(fit[[3]])
	}
	
	info$AUC[i] = AUC
}

info$Target = info$sample_name 
KO_info = info

load(tmpinf)
all_target = c("HepG2","0","All",mean(1-fit[[3]]),"All")

info = KO_info
info = apply(info,2, function(x) as.vector(x))
info = rbind(all_target, info)
info = as.data.frame(info)

row.names(info) = info$sample_name
info$Order_type = as.numeric(as.vector(info$Order_type))
info$AUC = as.numeric(as.vector(info$AUC))

info = info[order(info$Order_type),]
#info$Name = paste0(info$sample_name,"_",info$Target)
info$Name = info$Target
info$Name = factor(info$Name, levels = info$Name)
info$AUC = round(info$AUC,2)
AUC_info = info

library(randomForest)
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_shared_TF_data.Rda"
load(myinf1)

label = HepG2_data$tag


tag = AUC_info$AUC < 0.77
target = AUC_info$Target[tag==0]
target = as.vector(target)
tag = grep("All",target)
target = target[-tag]

tag = which(colnames(HepG2_data) %in% target)
HepG2_data = HepG2_data[,-tag]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_shared_core_TF_data.Rda"
save(HepG2_data, file = myoutf)

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(HepG2_data) 
HepG2_fit = fit

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_shared_core_TF_fit.Rda"
save(HepG2_fit,file = myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_shared_core_TF_model.Rda"
HepG2_model <- randomForest(tag~., data=HepG2_data, ntree=500,sampsize=c(2000,2000),importance=T)
save(HepG2_model,file = myoutf)


[3]GM12878 model in HepG2
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_shared_TF_data.Rda"
load(myinf1)

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_core_TF_model.Rda"
load(myinf2)
#########################
#GM12878 model in HepG2
########################
res1 = HepG2_data

#take the shared TF
myinf3 =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_core_TF_data.Rda"
load(myinf3)

com = intersect(colnames(GM12878_data),colnames(HepG2_data))
HepG2_data = HepG2_data[,com]

res1 = HepG2_data
HepG2_pred = predict(GM12878_model, res1[,-1], type="prob")
tmp = data.frame(res1[,1], HepG2_pred)
res = tmp
	
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
GM12878_model_in_HepG2_auc = 1- myauc

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_model_in_HepG2_core_TF_AUC.Rda"
save(GM12878_model_in_HepG2_auc, file=myoutf)

[4]HepG2 model in GM12878
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_TF_data.Rda"
load(myinf1)

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_shared_core_TF_model.Rda"
load(myinf2)
#########################
#HepG2 model in GM12878
########################
res1 = HepG2_data

#take the shared TF
myinf3 =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_shared_core_TF_data.Rda"
load(myinf3)

com = intersect(colnames(GM12878_data),colnames(HepG2_data))
GM12878_data = GM12878_data[,com]

res1 = GM12878_data
GM12878_pred = predict(HepG2_model, res1[,-1], type="prob")
tmp = data.frame(res1[,1], GM12878_pred)
res = tmp
	
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
HepG2_model_in_GM12878_auc = 1- myauc

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_model_in_GM12878_core_TF_AUC.Rda"
save(GM12878_model_in_HepG2_auc, file=myoutf)
