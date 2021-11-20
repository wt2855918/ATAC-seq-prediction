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

library(randomForest)
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/HepG2_Histone_integration_AUC.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_Histone_integration_AUC.Rda"
myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/mESC/mESC_Histone_integration_AUC.Rda"
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
#######################
load(myinf3)
mESC_info = info
mESC_data = data

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
mESC_info = final_info

row.names(GM12878_info) = GM12878_info$Experiment.target
row.names(HepG2_info) = HepG2_info$Experiment.target
row.names(mESC_info) = mESC_info$Experiment.target


res = matrix(0,11,3)
res = as.data.frame(res)
row.names(res) = row.names(GM12878_info)
colnames(res) = c("GM12878","HepG2","mESC")

res[row.names(GM12878_info),"GM12878"] = GM12878_info$AUC
res[row.names(HepG2_info),"HepG2"] = HepG2_info$AUC
res[row.names(mESC_info),"mESC"] = mESC_info$AUC

avg = apply(res[,1:2],1,mean)
avg = avg[order(avg,decreasing=T)]

res[names(avg),]


          GM12878     HepG2      mESC
H3K4me2  0.7057788 0.7527850 0.0000000
H2AFZ    0.7535010 0.6867868 0.0000000
H3K4me3  0.7035805 0.7268595 0.6068702
H3K36me3 0.5772770 0.8127395 0.5890053
H3K9ac   0.6762858 0.6922295 0.0000000
H3K27me3 0.6462940 0.6687703 0.5467405
H3K27ac  0.6657615 0.6156523 0.0000000
H3K4me1  0.6328325 0.6079845 0.5245185
H4K20me1 0.6142467 0.5919868 0.0000000
H3K9me3  0.5841905 0.5413077 0.5724885
H3K79me2 0.5373890 0.5294535 0.0000000




















