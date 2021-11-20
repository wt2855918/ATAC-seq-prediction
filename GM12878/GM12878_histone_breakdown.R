[1]Analysis 
rm(list=ls())
library(randomForest)
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_integration_prediction.Rda")
com = intersect(colnames(data),row.names(final_info))
data = data[,com]
final_info = final_info[com,]
colnames(data) = final_info$Experiment.target
tag = c(rep(1,2000),rep(0,2000))
data = cbind(tag,data)
data$tag = factor(data$tag,levels =c("0","1"))

#relative importance
res = fit[[6]]
com = intersect(row.names(res),row.names(final_info))
res = res[com,]
final_info = final_info[com,]
row.names(res) = final_info$Experiment.target

res = apply(res,1,mean)
res = res[order(res,decreasing=T)]
tar = names(res)

for(i in 1 : (length(tar)-1))
{
	cat("\r",i)
	
	sam_for_taken = tar[1:i]
	tag = which(colnames(data) %in% sam_for_taken)
	test = data[,-tag]
	
	test  = test[complete.cases(test),]
	test[,"tag"] = factor(test[,"tag"], levels =c(0,1))

	if(i== (length(tar)-1))
	{
		source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_Uni_variable_10_CV.R")
		fit = ROC_uni(test)
		myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_KO_RF/GM12878_Order_",i,"_",tar[i],"_KO_histone.Rda")
		save(fit, file = myoutf)
		
	}else{
	
		source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")	
		fit = ROC(test) 
		rank = fit[[6]]
		rank = apply(rank,1,mean)
		rank = rank[order(rank,decreasing=T)]
		tar[i+1] = names(rank)[1]	
		
		myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_KO_RF/GM12878_Order_",i,"_",tar[i],"_KO_histone.Rda")
		save(fit, file = myoutf)
		}
	
	
	#begin to calculate the rank
}

[2]Begin to plot
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_KO_RF/"
files = list.files(mydir)

tmpinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_integration_prediction.Rda"

my_nam = gsub("_KO_histone.Rda","",files)

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
c1 <- brewer.pal(10,"RdBu")[3]

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_2/Histone_KO_each.pdf"
pdf(myoutf1, width= 4.5, height= 3)
p <- position_dodge(0.1)
p1 <- ggplot(data=info, aes(x=Name, y=AUC, fill=c1)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_fill_manual(values=c(c1))+guides(fill=FALSE)+scale_y_continuous(expand = c(0,0),limits = c(0, 0.87)) + geom_text(aes(label=AUC), vjust=1.6, color="white", size=2)
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(angle=45,size= 7.25,hjust=1),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank()) + ylab("AUC score") + geom_hline(yintercept=0.5, linetype="dashed", color = "red", size=0.5)
p10 <- p9 +ggtitle(" ") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
p10 <- p10 + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
p10
dev.off()









