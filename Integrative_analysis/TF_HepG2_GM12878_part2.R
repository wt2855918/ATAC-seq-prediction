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

[1]Use GM12878 Specific to predict HepG2

library(randomForest)
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_integration_prediction.Rda"
load(myinf1)
HepG2_info = final_info
HepG2_data = data

com = intersect(colnames(HepG2_data),row.names(HepG2_info))
HepG2_data = HepG2_data[,com]
HepG2_info = HepG2_info[com,]
colnames(HepG2_data) = HepG2_info$Experiment.target

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda"
load(myinf2)
GM12878_info = final_info
GM12878_data = data
com = intersect(colnames(GM12878_data),row.names(GM12878_info))
GM12878_data = GM12878_data[,com]
GM12878_info = GM12878_info[com,]
colnames(GM12878_data) = GM12878_info$Experiment.target

com = intersect(colnames(GM12878_data),colnames(HepG2_data))
GM12878_data = GM12878_data[,com]
HepG2_data = HepG2_data[,com]

row.names(GM12878_info) = GM12878_info$Experiment.target
row.names(HepG2_info) = HepG2_info$Experiment.target

com = intersect(row.names(GM12878_info),row.names(HepG2_info))
GM12878_info = GM12878_info[com,]
HepG2_info = HepG2_info[com,]

#load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_sp_TF.Rda")
tar3 = c("TRIM22", "IKZF1", "ATF2")

GM12878_data = GM12878_data[,tar3]
HepG2_data = HepG2_data[,tar3]

GM12878_info = GM12878_info[tar3,]
HepG2_info = HepG2_info[tar3,]

row.names(GM12878_info) = GM12878_info$File.accession
row.names(HepG2_info) = HepG2_info$File.accession

###########################
#GM12878 self preidction
##############################
tag = c(rep(1,2000),rep(0,2000))
GM12878_data = cbind(tag,GM12878_data)
GM12878_data$tag = factor(GM12878_data$tag, levels = c("0","1"))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(GM12878_data) 
GM12878_fit = fit

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_sp_TF_fit.Rda"
save(fit,file = myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_sp_TF_data.Rda"
save(GM12878_data,file = myoutf)

#save the model for further usage
GM12878_model <- randomForest(tag~., data=GM12878_data, ntree=500,sampsize=c(2000,2000),importance=T)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_sp_TF_Model.Rda"
save(GM12878_model,file = myoutf)

#########################
#GM12878 model in HepG2
########################
res1 = HepG2_data
tag = c(rep(1,2000),rep(0,2000))
res1 = cbind(tag,res1)
res1 = as.data.frame(res1)
res1$tag = factor(res1$tag, levels =c(0,1))

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

GM12878_AUC = mean(1-fit[[3]])
HepG2_AUC = GM12878_model_in_HepG2_auc = 0.5


data = rbind(GM12878_AUC, HepG2_AUC)
data = as.data.frame(data)
data$target = c("GM12878","HepG2")
data$V1 = round(data$V1,2)
data$target= factor(data$target, levels = c("GM12878","HepG2") )
c1 <- brewer.pal(10,"RdBu")[3]

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_5/GM12878_sp_TF_HepG2.pdf"
pdf(myoutf1, width= 3, height= 3)
p <- position_dodge(0.1)
p1 <- ggplot(data=data, aes(x=target, y=V1, fill=c1)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5) + scale_fill_manual(values=c(c1))+guides(fill=FALSE)+scale_y_continuous(expand = c(0,0),limits = c(0, 0.85)) + geom_text(aes(label=V1), vjust=1.6, color="white", size=3.5)
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(angle=45,size= 7.25),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank()) + ylab("AUC score") + geom_hline(yintercept=0.5, linetype="dashed", color = "red", size=0.5)
p10 <- p9 +ggtitle(" ") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
p10 <- p10 + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
p10
dev.off()

[2] HepG2 Specific to predict GM12878
rm(list=ls())
library(randomForest)
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_integration_prediction.Rda"
load(myinf1)
HepG2_info = final_info
HepG2_data = data

com = intersect(colnames(HepG2_data),row.names(HepG2_info))
HepG2_data = HepG2_data[,com]
HepG2_info = HepG2_info[com,]
colnames(HepG2_data) = HepG2_info$Experiment.target

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda"
load(myinf2)
GM12878_info = final_info
GM12878_data = data
com = intersect(colnames(GM12878_data),row.names(GM12878_info))
GM12878_data = GM12878_data[,com]
GM12878_info = GM12878_info[com,]
colnames(GM12878_data) = GM12878_info$Experiment.target

com = intersect(colnames(GM12878_data),colnames(HepG2_data))
GM12878_data = GM12878_data[,com]
HepG2_data = HepG2_data[,com]

row.names(GM12878_info) = GM12878_info$Experiment.target
row.names(HepG2_info) = HepG2_info$Experiment.target

com = intersect(row.names(GM12878_info),row.names(HepG2_info))
GM12878_info = GM12878_info[com,]
HepG2_info = HepG2_info[com,]

raw.GM12878_data = GM12878_data
raw.HepG2_data = HepG2_data

GM12878_data = raw.GM12878_data
HepG2_data = raw.HepG2_data

#load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_sp_TF.Rda")
tar1 = c("REST", "USF1", "CEBPB")
GM12878_data = GM12878_data[,tar1]
HepG2_data = HepG2_data[,tar1]

GM12878_info = GM12878_info[tar1,]
HepG2_info = HepG2_info[tar1,]

row.names(GM12878_info) = GM12878_info$File.accession
row.names(HepG2_info) = HepG2_info$File.accession

###########################
#HepG2 self preidction
##############################
tag = c(rep(1,2000),rep(0,2000))
HepG2_data = cbind(tag,HepG2_data)
HepG2_data$tag = factor(HepG2_data$tag, levels = c("0","1"))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(HepG2_data) 
HepG2_fit = fit
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/HepG2_sp_TF_fit.Rda"
save(fit,file = myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_sp_TF_data.Rda"
save(HepG2_data,file = myoutf)

#save the model for further usage
HepG2_model <- randomForest(tag~., data=HepG2_data, ntree=500,sampsize=c(2000,2000),importance=T)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_sp_TF_Model.Rda"
save(HepG2_model,file = myoutf)

#########################
#HepG2 model in GM12878
########################
res2 = GM12878_data
tag = c(rep(1,2000),rep(0,2000))
res2 = cbind(tag,res2)
res2 = as.data.frame(res2)
res2$tag = factor(res2$tag, levels =c(0,1))

GM12878_pred = predict(HepG2_model, res2[,-1], type="prob")
tmp = data.frame(res2[,1], GM12878_pred)
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

HepG2_AUC = mean(1-fit[[3]])
GM12878_AUC = HepG2_model_in_GM12878_auc

data = rbind(GM12878_AUC, HepG2_AUC)
data = as.data.frame(data)
data$target = c("GM12878","HepG2")
data$V1 = round(data$V1,2)
data$target= factor(data$target, levels = c("HepG2","GM12878") )
c1 <- brewer.pal(10,"RdBu")[8]

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_5/HepG2_sp_TF_HepG2.pdf"
pdf(myoutf1, width= 3, height= 3)
p <- position_dodge(0.1)
p1 <- ggplot(data=data, aes(x=target, y=V1, fill=c1)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5) + scale_fill_manual(values=c(c1))+guides(fill=FALSE)+scale_y_continuous(expand = c(0,0),limits = c(0, 0.85)) + geom_text(aes(label=V1), vjust=1.6, color="white", size=3.5)
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(angle=45,size= 7.25),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank()) + ylab("AUC score") + geom_hline(yintercept=0.5, linetype="dashed", color = "red", size=0.5)
p10 <- p9 +ggtitle(" ") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
p10 <- p10 + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
p10
dev.off()

[3]GM12878 AND HepG2
library(randomForest)
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_integration_prediction.Rda"
load(myinf1)
HepG2_info = final_info
HepG2_data = data

com = intersect(colnames(HepG2_data),row.names(HepG2_info))
HepG2_data = HepG2_data[,com]
HepG2_info = HepG2_info[com,]
colnames(HepG2_data) = HepG2_info$Experiment.target

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda"
load(myinf2)
GM12878_info = final_info
GM12878_data = data
com = intersect(colnames(GM12878_data),row.names(GM12878_info))
GM12878_data = GM12878_data[,com]
GM12878_info = GM12878_info[com,]
colnames(GM12878_data) = GM12878_info$Experiment.target

com = intersect(colnames(GM12878_data),colnames(HepG2_data))
GM12878_data = GM12878_data[,com]
HepG2_data = HepG2_data[,com]

row.names(GM12878_info) = GM12878_info$Experiment.target
row.names(HepG2_info) = HepG2_info$Experiment.target

com = intersect(row.names(GM12878_info),row.names(HepG2_info))
GM12878_info = GM12878_info[com,]
HepG2_info = HepG2_info[com,]


#load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_sp_TF.Rda")
tar2 = c("CREM","BHLHE40","YY1","ARID3A")
GM12878_data = GM12878_data[,tar2]
HepG2_data = HepG2_data[,tar2]

GM12878_info = GM12878_info[tar2,]
HepG2_info = HepG2_info[tar2,]

row.names(GM12878_info) = GM12878_info$File.accession
row.names(HepG2_info) = HepG2_info$File.accession

###########################
#GM12878 self preidction
##############################
tag = c(rep(1,2000),rep(0,2000))
GM12878_data = cbind(tag,GM12878_data)
GM12878_data$tag = factor(GM12878_data$tag, levels = c("0","1"))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(GM12878_data) 
GM12878_fit = fit

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_core_shared_TF_fit.Rda"
save(fit,file = myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_core_shared_TF_data.Rda"
save(GM12878_data,file = myoutf)

#save the model for further usage
GM12878_model <- randomForest(tag~., data=GM12878_data, ntree=500,sampsize=c(2000,2000),importance=T)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_core_shared_TF_Model.Rda"
save(GM12878_model,file = myoutf)


###########################
#HepG2 self preidction
##############################
tag = c(rep(1,2000),rep(0,2000))
HepG2_data = cbind(tag,HepG2_data)
HepG2_data$tag = factor(HepG2_data$tag, levels = c("0","1"))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(HepG2_data) 
HepG2_fit = fit
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/HepG2_core_shared_TF_fit.Rda"
save(fit,file = myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_core_shared_TF_data.Rda"
save(HepG2_data,file = myoutf)

#save the model for further usage
HepG2_model <- randomForest(tag~., data=HepG2_data, ntree=500,sampsize=c(2000,2000),importance=T)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_core_shared_TF_Model.Rda"
save(HepG2_model,file = myoutf)


res1 = HepG2_data[,2:ncol(HepG2_data)]
res2 = GM12878_data[,2:ncol(GM12878_data)]
#########################
#HepG2 model in GM12878
########################
tag = c(rep(1,2000),rep(0,2000))
res2 = cbind(tag,res2)
res2 = as.data.frame(res2)
res2$tag = factor(res2$tag, levels =c(0,1))

GM12878_pred = predict(HepG2_model, res2[,-1], type="prob")
tmp = data.frame(res2[,1], GM12878_pred)
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

#########################
#GM12878 model in HepG2
########################
tag = c(rep(1,2000),rep(0,2000))
res1 = cbind(tag,res1)
res1 = as.data.frame(res1)
res1$tag = factor(res1$tag, levels =c(0,1))

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

#HepG2 figure
HepG2_AUC = mean(1-HepG2_fit[[3]])

data = rbind(GM12878_model_in_HepG2_auc, mean(HepG2_AUC))
data = as.data.frame(data)
data$target = c("GM12878_Model","HepG2")
data$V1 = round(data$V1,2)
data$target= factor(data$target, levels = c("GM12878_Model","HepG2") )
c1 <- brewer.pal(10,"RdBu")[3]

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_5/GM12878_core_TF_Model_in_HepG2_AUC.pdf"
pdf(myoutf1, width= 3, height= 3)
p <- position_dodge(0.1)
p1 <- ggplot(data=data, aes(x=target, y=V1, fill=c1)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5) + scale_fill_manual(values=c(c1))+guides(fill=FALSE)+scale_y_continuous(expand = c(0,0),limits = c(0, 0.85)) + geom_text(aes(label=V1), vjust=1.6, color="white", size=3.5)
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(angle=45,size= 7.25),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank()) + ylab("AUC score") + geom_hline(yintercept=0.5, linetype="dashed", color = "red", size=0.5)
p10 <- p9 +ggtitle(" ") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
p10 <- p10 + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
p10
dev.off()

#GM12878 figure
GM12878_AUC = mean(1-GM12878_fit[[3]])

data = rbind( HepG2_model_in_GM12878_auc, mean(GM12878_AUC))
data = as.data.frame(data)
data$target = c("HepG2_Model","GM12878")
data$V1 = round(data$V1,2)
data$target= factor(data$target, levels = c("HepG2_Model","GM12878") )
c1 <- brewer.pal(10,"RdBu")[8]

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_5/HepG2_core_TF_Model_in_GM12878_AUC.pdf"
pdf(myoutf1, width= 3, height= 3)
p <- position_dodge(0.1)
p1 <- ggplot(data=data, aes(x=target, y=V1, fill=c1)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5) + scale_fill_manual(values=c(c1))+guides(fill=FALSE)+scale_y_continuous(expand = c(0,0),limits = c(0, 0.85)) + geom_text(aes(label=V1), vjust=1.6, color="white", size=3.5)
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(angle=45,size= 7.25),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank()) + ylab("AUC score") + geom_hline(yintercept=0.5, linetype="dashed", color = "red", size=0.5)
p10 <- p9 +ggtitle(" ") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
p10 <- p10 + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
p10
dev.off()