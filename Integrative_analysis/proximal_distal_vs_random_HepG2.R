[1]HepG2 model TF in proximal/distal vs random
rm(list=ls())
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_integration_prediction.Rda")
res1 = data

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_TF_integration_Core.Rda")

com = intersect(colnames(res1),colnames(res2))
res2 = res2[,com]

data = res2
tag = data$tag == "0"
data = data[tag,]
data = data[,2:ncol(data)]

#change to histone names
#change the name 
com = intersect(colnames(data), row.names(final_info))
data = data[,com]
final_info = final_info[com,]
colnames(data) = final_info$Experiment.target
tag = sample(nrow(data), 2000)
data = data[tag,]
tag = rep(0,2000)
data = cbind(tag,data)
data_neg = data

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_prediction_promoter_distal.Rda")

###############
#proximal data
################
Annotated_data = data
proximal_data_pos = Annotated_data[Annotated_data$tag==1,]
distal_data_pos = Annotated_data[Annotated_data$tag==0,]
proximal_data = rbind(proximal_data_pos, data_neg)

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_prediction_model.Rda")

res = NULL
tmp = predict(model, proximal_data[,-1], type="prob")
tmp = data.frame(proximal_data[,1], tmp)
res = rbind(res, tmp)

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

Model_proximal_xx = xx
Model_proximal_yy = yy
Model_proximal_AUC = 1- myauc

###############
#distal data
################
distal_data = rbind(distal_data_pos, data_neg)
distal_data$tag = c(rep(1,2000),rep(0,2000))
distal_data$tag = factor(distal_data$tag, levels=c(0,1))

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_prediction_model.Rda")

res = NULL
tmp = predict(model, distal_data[,-1], type="prob")
tmp = data.frame(distal_data[,1], tmp)
res = rbind(res, tmp)

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

Model_distal_xx = xx
Model_distal_yy = yy
Model_distal_AUC = 1- myauc

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_prediction_promoter_distal_vs_random.Rda"
save.image(file= myoutf)

###############
#merge data
################
merge_data = rbind(proximal_data_pos, distal_data_pos, data_neg)
merge_data$tag = c(rep(1,4000),rep(0,2000))
merge_data$tag = factor(merge_data$tag, levels=c(0,1))

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_prediction_model.Rda")

res = NULL
tmp = predict(model, merge_data[,-1], type="prob")
tmp = data.frame(merge_data[,1], tmp)
res = rbind(res, tmp)

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

Model_distal_xx = xx
Model_distal_yy = yy
Model_distal_AUC = 1- myauc


[2]GM12878 model HM in proximal/distal vs random
rm(list=ls())
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_integration_prediction.Rda")
res2 = data

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_Histone_integration_Core.Rda")

com = intersect(colnames(res1),colnames(res2))
res1 = res1[,com]

data = res1
tag = data$tag == "0"
data = data[tag,]
data = data[,2:ncol(data)]

#change to histone names
#change the name 
com = intersect(colnames(data), row.names(final_info))
data = data[,com]
final_info = final_info[com,]
colnames(data) = final_info$Experiment.target
tag = sample(nrow(data), 2000)
data = data[tag,]
tag = rep(0,2000)
data = cbind(tag,data)
data_neg = data

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_prediction_promoter_distal.Rda")

###############
#proximal data
################
Annotated_data = data
proximal_data_pos = Annotated_data[Annotated_data$tag==1,]
distal_data_pos = Annotated_data[Annotated_data$tag==0,]
proximal_data = rbind(proximal_data_pos, data_neg)

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_prediction_model.Rda")

res = NULL
tmp = predict(model, proximal_data[,-1], type="prob")
tmp = data.frame(proximal_data[,1], tmp)
res = rbind(res, tmp)

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

Model_proximal_xx = xx
Model_proximal_yy = yy
Model_proximal_AUC = 1- myauc

###############
#distal data
################
distal_data = rbind(distal_data_pos, data_neg)
distal_data$tag = c(rep(1,2000),rep(0,2000))
distal_data$tag = factor(distal_data$tag, levels=c(0,1))

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_prediction_model.Rda")

res = NULL
tmp = predict(model, distal_data[,-1], type="prob")
tmp = data.frame(distal_data[,1], tmp)
res = rbind(res, tmp)

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

Model_distal_xx = xx
Model_distal_yy = yy
Model_distal_AUC = 1- myauc

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_prediction_promoter_distal_vs_random.Rda"
save.image(file= myoutf)

#######################################
#################################
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_prediction_promoter_distal_vs_random.Rda"
load(myinf1)
TF_auc = Model_proximal_AUC

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_prediction_promoter_distal_vs_random.Rda"
load(myinf2)
histone_auc = Model_proximal_AUC

data = rbind(TF_auc, histone_auc)

data = as.data.frame(data)
data$target = c("TF","Histone")
colnames(data)[1] = c("AUC")
data$AUC = round(data$AUC,2)

#c1 <- brewer.pal(10,"RdBu")[3]

c1 <- brewer.pal(10,"Set1")[1]
c2 <- brewer.pal(10,"Set1")[2]
color <- c(c1,c2)
data$target = factor(data$target,levels=c("TF","Histone"))

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/HepG2/TF_histone_other_Model_proximal_vs_random_prediction_barplot.pdf"
pdf(myoutf1, width= 3, height= 3)
p <- position_dodge(0.1)
p1 <- ggplot(data=data, aes(x=target, y=AUC, fill=color)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5) + scale_fill_manual(values=c(c2,c1))+guides(fill=FALSE)+scale_y_continuous(expand = c(0,0),limits = c(0, 1)) + geom_text(aes(label=AUC), vjust=1.6, color="white", size=3.5)
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(angle=45,size= 7.25),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank()) + ylab("AUC score") + geom_hline(yintercept=0.5, linetype="dashed", color = "red", size=0.5)
p10 <- p9 +ggtitle(" ") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
p10 <- p10 + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
p10
dev.off()


#################################
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_prediction_promoter_distal_vs_random.Rda"
load(myinf1)
TF_auc = Model_distal_AUC

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_prediction_promoter_distal_vs_random.Rda"
load(myinf2)
histone_auc = Model_distal_AUC

data = rbind(TF_auc, histone_auc)

data = as.data.frame(data)
data$target = c("TF","Histone")
colnames(data)[1] = c("AUC")
data$AUC = round(data$AUC,2)

#c1 <- brewer.pal(10,"RdBu")[3]

c1 <- brewer.pal(10,"Set1")[1]
c2 <- brewer.pal(10,"Set1")[2]
color <- c(c1,c2)
data$target = factor(data$target,levels=c("TF","Histone"))

#myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/HepG2/TF_histone_other_Model_distal_vs_random_prediction_barplot.pdf"
myoutf1 = "/ihome/yanding/HepG2_TF_histone_other_Model_distal_vs_random_prediction_barplot.pdf"
pdf(myoutf1, width= 3, height= 3)
p <- position_dodge(0.1)
p1 <- ggplot(data=data, aes(x=target, y=AUC, fill=color)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5) + scale_fill_manual(values=c(c2,c1))+guides(fill=FALSE)+scale_y_continuous(expand = c(0,0),limits = c(0, 1)) + geom_text(aes(label=AUC), vjust=1.6, color="white", size=3.5)
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(angle=45,size= 7.25),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank()) + ylab("AUC score") + geom_hline(yintercept=0.5, linetype="dashed", color = "red", size=0.5)
p10 <- p9 +ggtitle(" ") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
p10 <- p10 + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
p10
dev.off()

