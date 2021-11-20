
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

[1]Use GM12878 to predict HepG2

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

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_TF_fit.Rda"
save(fit,file = myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_TF_data.Rda"
save(GM12878_data,file = myoutf)

#save the model for further usage
GM12878_model <- randomForest(tag~., data=GM12878_data, ntree=500,sampsize=c(2000,2000),importance=T)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_TF_Model.Rda"
save(GM12878_model,file = myoutf)

#begin to plot the GM12878 figure relative importance
importance = fit[[6]]
rank = apply(importance,1,mean)
rank = rank[order(rank,decreasing=T)]

for(i in 1 : ncol(importance))
{
	cat("\r",i)
	
	xx = importance[,i]
	xx = as.data.frame(xx)
	xx$target = row.names(xx)
	
	if(i == 1)
	{
		res = xx
	
	}else{
		
		res = rbind(xx,res)
	
	}

}
data = res
colnames(data) <- c("EDS", "group")
data$EDS <- as.numeric(as.vector(data$EDS))

library(RColorBrewer)
#c1 <- brewer.pal(10,"RdBu")[3]
c2 <- brewer.pal(10,"RdBu")[3]

tar <- c("EDS")
plot_gene <- function(data,tar)
{
 
 
  data[,"group"] <- factor(data[,"group"],levels=names(rank))
  p1 <- ggplot(data, aes(group,data[,tar],fill=data[,"group"]))+
    geom_boxplot(outlier.shape = NA,position="dodge",notch = F,width=0.5,aes(fill=factor(group)))
  p2 <- p1 +ylab("Importance")
  p3 <- p2 + theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p3
}


myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_5/GM12878_shared_TF_importance_bxplot.pdf"
pdf(myoutf1, width= 4.5, height= 3)
p <- position_dodge(0.1)
p1 <- plot_gene(data,tar) + scale_fill_manual(values=rep(c2,36))+guides(fill=FALSE) + scale_y_continuous(breaks= pretty(data$EDS,n=4))
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(size= 7.25,angle=90),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank())
p10 <- p9 +ggtitle("GM12878") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
#p10 <- p10+ annotate("text",x=1,y=py[length(py)-2],label=paste0("p=",pval),size=2.5,)
p10
dev.off()

###########################
#HepG2 self preidction
##############################
tag = c(rep(1,2000),rep(0,2000))
HepG2_data = cbind(tag,HepG2_data)
HepG2_data$tag = factor(HepG2_data$tag, levels = c("0","1"))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(HepG2_data) 
HepG2_fit = fit
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/HepG2_shared_TF_fit.Rda"
save(fit,file = myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_shared_TF_data.Rda"
save(HepG2_data,file = myoutf)

#save the model for further usage
HepG2_model <- randomForest(tag~., data=HepG2_data, ntree=500,sampsize=c(2000,2000),importance=T)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_shared_TF_Model.Rda"
save(HepG2_model,file = myoutf)

#begin to plot the GM12878 figure relative importance
importance = fit[[6]]
rank = apply(importance,1,mean)
rank = rank[order(rank,decreasing=T)]

for(i in 1 : ncol(importance))
{
	cat("\r",i)
	
	xx = importance[,i]
	xx = as.data.frame(xx)
	xx$target = row.names(xx)
	
	if(i == 1)
	{
		res = xx
	
	}else{
		
		res = rbind(xx,res)
	
	}

}
data = res
colnames(data) <- c("EDS", "group")
data$EDS <- as.numeric(as.vector(data$EDS))

library(RColorBrewer)
#c1 <- brewer.pal(10,"RdBu")[3]
c2 <- brewer.pal(10,"RdBu")[8]

tar <- c("EDS")
plot_gene <- function(data,tar)
{
 
 
  data[,"group"] <- factor(data[,"group"],levels=names(rank))
  p1 <- ggplot(data, aes(group,data[,tar],fill=data[,"group"]))+
    geom_boxplot(outlier.shape = NA,position="dodge",notch = F,width=0.5,aes(fill=factor(group)))
  p2 <- p1 +ylab("Importance")
  p3 <- p2 + theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p3
}


myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_5/HepG2_shared_TF_importance_bxplot.pdf"
pdf(myoutf1, width= 4.5, height= 3)
p <- position_dodge(0.1)
p1 <- plot_gene(data,tar) + scale_fill_manual(values=rep(c2,36))+guides(fill=FALSE) + scale_y_continuous(breaks= pretty(data$EDS,n=4))
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(size= 7.25,angle=90),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank())
p10 <- p9 +ggtitle("HepG2") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
#p10 <- p10+ annotate("text",x=1,y=py[length(py)-2],label=paste0("p=",pval),size=2.5,)
p10
dev.off()

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

#############################################
#GM12878 signal in HepG2 using GM12878 model
#############################################
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/"
files = paste0(row.names(GM12878_info),".tab")

res = matrix(0,4000,36)
row.names(res) = row.names(HepG2_data)
colnames(res) = GM12878_info$Experiment.target
res = as.data.frame(res)

for(i in 1 : length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir,files[i])
	tmp = read.table(tmpinf,sep="\t",quote=NULL)
	tmp$V3 = tmp$V3 - 1
	row.names(tmp) = paste0(tmp$V1,"_",tmp$V2,"_",tmp$V3)
	res[,i] = tmp[row.names(HepG2_data),"V4"]
	
}
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_TF_signal_HepG2_validation_region.Rda"
save(res,file =myoutf)

tag = c(rep(1,2000),rep(0,2000))
res = cbind(tag,res)
res$tag = factor(res$tag, levels =c(0,1))
res = res[complete.cases(res),]

HepG2_pred_GM12878_sig = predict(GM12878_model, res[,-1], type="prob")
tmp = data.frame(res[,1], HepG2_pred_GM12878_sig)
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
GM12878_model_in_HepG2_GM12878_sig_auc = 1- myauc

#############################################
#HepG2 signal in GM12878 using HepG2 model
#############################################
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
files = paste0(row.names(HepG2_info),".tab")

res = matrix(0,4000,36)
row.names(res) = row.names(GM12878_data)
colnames(res) = HepG2_info$target
res = as.data.frame(res)

for(i in 1 : length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir,files[i])
	tmp = read.table(tmpinf,sep="\t",quote=NULL)
	tmp$V3 = tmp$V3 - 1
	row.names(tmp) = paste0(tmp$V1,"_",tmp$V2,"_",tmp$V3)
	res[,i] = tmp[row.names(GM12878_data),"V4"]
	
}
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_TF_signal_GM12878_validation_region.Rda"
save(res,file =myoutf)

tag = c(rep(1,2000),rep(0,2000))
res = cbind(tag,res)
res$tag = factor(res$tag, levels =c(0,1))
res = res[complete.cases(res),]

HepG2_pred_GM12878_sig = predict(HepG2_model, res[,-1], type="prob")
tmp = data.frame(res[,1], HepG2_pred_GM12878_sig)
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
HepG2_model_in_GM12878_HepG2_sig_auc  = 1- myauc


#HepG2 figure
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/HepG2_shared_TF_fit.Rda")
HepG2_AUC = mean(1-fit[[3]])

data = rbind(GM12878_model_in_HepG2_GM12878_sig_auc, GM12878_model_in_HepG2_auc, mean(HepG2_AUC))
data = as.data.frame(data)
data$target = c("GM12878_Signal","GM12878_Model","HepG2")
data$V1 = round(data$V1,2)
data$target= factor(data$target, levels = c("GM12878_Signal","GM12878_Model","HepG2") )
c1 <- brewer.pal(10,"RdBu")[3]

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_5/GM12878_TF_Model_in_HepG2_AUC.pdf"
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
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_TF_fit.Rda")
GM12878_AUC = mean(1-fit[[3]])

data = rbind(HepG2_model_in_GM12878_HepG2_sig_auc, HepG2_model_in_GM12878_auc, mean(GM12878_AUC))
data = as.data.frame(data)
data$target = c("HepG2_Signal","HepG2_Model","GM12878")
data$V1 = round(data$V1,2)
data$target= factor(data$target, levels = c("HepG2_Signal","HepG2_Model","GM12878") )
c1 <- brewer.pal(10,"RdBu")[8]

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_5/HepG2_Histone_Model_in_GM12878_AUC.pdf"
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


#save the information
library(randomForest)
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_integration_prediction.Rda"
load(myinf1)
HepG2_info = final_info
HepG2_data = data

com = intersect(colnames(HepG2_data),row.names(HepG2_info))
HepG2_data = HepG2_data[,com]
HepG2_info = HepG2_info[com,]
colnames(HepG2_data) = HepG2_info$Experiment.target

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/HepG2_TF_each_AUC.xls"
write.table(HepG2_info,myoutf,sep="\t",quote=F)


myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda"
load(myinf2)
GM12878_info = final_info
GM12878_data = data

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/GM12878_TF_each_AUC.xls"
write.table(GM12878_info,myoutf,sep="\t",quote=F)

row.names(HepG2_info) = HepG2_info$Experiment.target
row.names(GM12878_info) = GM12878_info$Experiment.target

tar = c("TRIM22","ELF1","IKZF1")
HepG2_info[tar,]
GM12878_info[tar,]

HepG2_info[tar,]
       File.accession Experiment.target       AUC
TRIM22    ENCFF090OBV            TRIM22 0.5657428
ELF1      ENCFF573GDY              ELF1 0.7167902
IKZF1     ENCFF284DOX             IKZF1 0.5362040
> GM12878_info[tar,]
       File.accession Experiment.target       AUC
TRIM22    ENCFF896KWL            TRIM22 0.7395010
ELF1      ENCFF937HON              ELF1 0.7092840
IKZF1     ENCFF405CHV             IKZF1 0.7767895
> 
############################################3

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_TF_fit.Rda")
#begin to plot the GM12878 figure relative importance
importance = fit[[6]]
rank = apply(importance,1,mean)
rank = rank[order(rank,decreasing=T)]

for(i in 1 : ncol(importance))
{
	cat("\r",i)
	
	xx = importance[,i]
	xx = as.data.frame(xx)
	xx$target = row.names(xx)
	
	if(i == 1)
	{
		res = xx
	
	}else{
		
		res = rbind(xx,res)
	
	}

}
data = res
colnames(data) <- c("EDS", "group")
data$EDS <- as.numeric(as.vector(data$EDS))
data$group = factor(data,levels = names(rank))
GM12878_plot_data = data

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/HepG2_shared_TF_fit.Rda")
#begin to plot the GM12878 figure relative importance
importance = fit[[6]]
rank = apply(importance,1,mean)
rank = rank[order(rank,decreasing=T)]

for(i in 1 : ncol(importance))
{
	cat("\r",i)
	
	xx = importance[,i]
	xx = as.data.frame(xx)
	xx$target = row.names(xx)
	
	if(i == 1)
	{
		res = xx
	
	}else{
		
		res = rbind(xx,res)
	
	}

}
data = res
colnames(data) <- c("EDS", "group")
data$EDS <- as.numeric(as.vector(data$EDS))

library(RColorBrewer)
#c1 <- brewer.pal(10,"RdBu")[3]
c2 <- brewer.pal(10,"RdBu")[8]

tar <- c("EDS")
plot_gene <- function(data,tar)
{
 
 
  data[,"group"] <- factor(data[,"group"],levels=levels(GM12878_plot_data$group))
  p1 <- ggplot(data, aes(group,data[,tar],fill=data[,"group"]))+
    geom_boxplot(outlier.shape = NA,position="dodge",notch = F,width=0.5,aes(fill=factor(group)))
  p2 <- p1 +ylab("Importance")
  p3 <- p2 + theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p3
}


myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_5/HepG2_shared_TF_importance_bxplot.pdf"
pdf(myoutf1, width= 4.5, height= 3)
p <- position_dodge(0.1)
p1 <- plot_gene(data,tar) + scale_fill_manual(values=rep(c2,36))+guides(fill=FALSE) + scale_y_continuous(breaks= pretty(data$EDS,n=4))
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(size= 7.25,angle=90),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank())
p10 <- p9 +ggtitle("HepG2") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
#p10 <- p10+ annotate("text",x=1,y=py[length(py)-2],label=paste0("p=",pval),size=2.5,)
p10
dev.off()

#############################
rm(list=ls())
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_TF_fit.Rda")
#begin to plot the GM12878 figure relative importance
importance = fit[[6]]
res1 = importance


load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/HepG2_shared_TF_fit.Rda")
#begin to plot the GM12878 figure relative importance
importance = fit[[6]]
res2 = importance

com = intersect(row.names(res1),row.names(res2))
res1 = res1[com,]
res2 = res2[com,]

xx = apply(res1,1,mean)
yy = apply(res2,1,mean)

data = cbind(x=xx,y=yy)
data = as.data.frame(data)
data$group = row.names(data)
df = data
colnames(df)[1:2] <- c("x","y")
px <- pretty(df$x)
py <- pretty(df$y)

#correlation coefficient
fit <- cor(df$x,df$y, method="s",use="complete.obs")
coef <- round(fit,2)

library(RColorBrewer)
#c1 <- brewer.pal(10,"RdBu")[3]
c2 <- brewer.pal(10,"Set1")[2]



commonTheme = list(labs(color="Density",fill="Density",
                        x="GM12878 importance",
                        y="HepG2 importance"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_5/GM12878_TF_importance_vs_HepG2_importance.pdf"
pdf(myoutf1, width= 3, height= 3)
p <- ggplot(df, aes(x=x,y=y,label = group)) +   theme(legend.position="none")  +  geom_point(size=1,color=c2) + geom_text_repel(size=2)
p1 <- p + geom_smooth(method=lm,linetype=1,colour="grey70",se=F,size=1) + ggtitle("")
p2 <- p1 + theme(panel.grid.minor = element_line(colour="green", size=0.5),panel.grid.major = element_line(colour="green", size=0.5)) + theme_bw()
p3 <- p2 + guides(alpha="none",fill=F) + commonTheme + theme(aspect.ratio=1,panel.border = element_rect(colour = "black",size=1),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=1),axis.text=element_text(size= 7.25),legend.position="none") 
p3 <- p3 +  scale_x_continuous(breaks=px, limits=range(px)) + scale_y_continuous(breaks=py, limits=range(py))
p3 <- p3 + annotate("text",x=px[length(px)-1],y=py[length(py)-1],label=paste0("rho=",coef),size=1.5,)
p3 <- p3 + geom_hline(yintercept=15, linetype="dashed", color = "red", size=0.5) + geom_vline(xintercept=15, linetype="dashed", color = "red", size=0.5)

print(p3)
dev.off()


#############################
rm(list=ls())
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_TF_fit.Rda")
#begin to plot the GM12878 figure relative importance
importance = fit[[6]]
res1 = importance


load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/HepG2_shared_TF_fit.Rda")
#begin to plot the GM12878 figure relative importance
importance = fit[[6]]
res2 = importance

com = intersect(row.names(res1),row.names(res2))
res1 = res1[com,]
res2 = res2[com,]
xx = apply(res1,1,mean)
yy = apply(res2,1,mean)

data = cbind(x=xx,y=yy)
data = as.data.frame(data)
data$group = row.names(data)

tag1 = data$x<10&data$y>15
tag2 = data$x>15&data$y>15
tag3 = data$x>15&data$y<15

tar1 = row.names(data)[tag1]
tar2 = row.names(data)[tag2]
tar3 = row.names(data)[tag3]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_sp_TF.Rda"
save(tar1,file = myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/Shared_H_G_TF.Rda"
save(tar2,file = myoutf)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_sp_TF.Rda"
save(tar3,file = myoutf)







myinf1= "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/HepG2_TF_each_AUC.xls"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/GM12878_TF_each_AUC.xls"

res1 = read.table(myinf1,sep="\t",quote=NULL)
res2 = read.table(myinf2,sep="\t",quote=NULL)

row.names(res1) = res1$Experiment.target
row.names(res2) = res2$Experiment.target

com = intersect(row.names(res1),row.names(res2))
res1 = res1[com,]
res2 = res2[com,]

res = cbind(res1$AUC,res2$AUC)
row.names(res)= com
colnames(res) = c("H_AUC","G_AUC")
res = as.data.frame(res)





