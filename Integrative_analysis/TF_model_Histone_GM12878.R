
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
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_Histone_integration_AUC.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_Histone_integration_AUC.Rda"

load(myinf1)
HepG2_info = info
HepG2_data = data

data = data[,2:ncol(data)]
pos = seq(1,2000,1)
neg = seq(2001,4000,1)
res = score_auc(data,pos,neg)

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

HepG2_info = final_info

#########################
load(myinf2)
GM12878_info = info
GM12878_data = data

data = data[,2:ncol(data)]
pos = seq(1,2000,1)
neg = seq(2001,4000,1)
res = score_auc(data,pos,neg)

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

GM12878_info = final_info

colnames(HepG2_info) = colnames(GM12878_info) = c("file_accession","target")


row.names(HepG2_info) = HepG2_info$target
row.names(GM12878_info) = GM12878_info$target

com = intersect(row.names(HepG2_info),row.names(GM12878_info))
HepG2_info = HepG2_info[com,]
GM12878_info = GM12878_info[com,]

row.names(HepG2_info) = HepG2_info$file_accession
row.names(GM12878_info) = GM12878_info$file_accession

HepG2_data = HepG2_data[,row.names(HepG2_info)]
GM12878_data = GM12878_data[,row.names(GM12878_info)]

colnames(HepG2_data) = HepG2_info$target
colnames(GM12878_data) = GM12878_info$target

res1 = HepG2_data
res2 = GM12878_data

res1_pos = res1[1:2000,]
res1_neg = res1[2001:4000,]

res2_pos = res2[1:2000,]
res2_neg = res2[2001:4000,]

tag_pos = sample(seq(1,2000),1000)
tag_neg = sample(seq(1,2000),1000)

res1_pos = res1_pos[tag_pos,]
res1_neg = res1_neg[tag_neg,]

res2_pos = res2_pos[tag_pos,]
res2_neg = res2_neg[tag_neg,]

tag = c(rep(1,1000),rep(0,1000))
res1 = rbind(res1_pos, res1_neg)
res2 = rbind(res2_pos, res2_neg)

#########################
#Model training in HepG2
########################
res1 = cbind(tag,res1)
res1 = as.data.frame(res1)
res1$tag = factor(res1$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(res1)

HepG2_AUC = 1- fit[[3]]

HepG2_model = randomForest(tag~., data=res1, ntree=500,importance=T)

#########################
#Model training in GM12878
########################
res2 = cbind(tag,res2)
res2 = as.data.frame(res2)
res2$tag = factor(res2$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(res2)

GM12878_AUC = 1- fit[[3]]

GM12878_model = randomForest(tag~., data=res2, ntree=500,importance=T)

#########################
#HepG2 model in GM12878
########################
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

data = rbind(GM12878_model_in_HepG2_auc, mean(HepG2_AUC))
data = as.data.frame(data)
data$target = c("GM12878_Model","HepG2")
data$V1 = round(data$V1,2)
c1 <- brewer.pal(10,"RdBu")[3]

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/GM12878_Histone_Model_in_HepG2_AUC.pdf"
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

data = rbind(HepG2_model_in_GM12878_auc, mean(GM12878_AUC))
data = as.data.frame(data)
data$target = c("HepG2_Model","GM12878")
data$V1 = round(data$V1,2)
c1 <- brewer.pal(10,"RdBu")[3]

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/HepG2/HepG2_Histone_Model_in_GM12878_AUC.pdf"
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







	