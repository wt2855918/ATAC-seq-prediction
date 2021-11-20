load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_TF_integration_Core.Rda")
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/HepG2_training_and_validation_pos.Rda"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/HepG2_training_and_validation_neg.Rda"

#define the genome ranges
peakDF2GRanges <- function(peak.df) {
    peak.gr=GRanges(seqnames=peak.df[,1],
        ranges=IRanges(peak.df[,2], peak.df[,3]))
    cn <- colnames(peak.df)
    if (length(cn) > 3) {
        for (i in 4:length(cn)) {
            mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
        }
    }
    return(peak.gr)
}



[1]Annotate the positive regions
rm(list=ls())
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(GenomicRanges)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_TF_integration_Core.Rda")
tag = res2$tag == "1"
res2 = res2[tag==1,]

chrom = as.vector(sapply(row.names(res2), function(x) strsplit(x,"_")[[1]][1]))
start = as.vector(sapply(row.names(res2), function(x) strsplit(x,"_")[[1]][2]))
end = as.vector(sapply(row.names(res2), function(x) strsplit(x,"_")[[1]][3]))

res = cbind(chrom,start,end)
res = as.data.frame(res)
res$start = as.numeric(as.vector(res$start))
res$end = as.numeric(as.vector(res$end))
res = peakDF2GRanges(res)
peakAnno <- annotatePeak(res, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")
info = peakAnno@anno
info = as.data.frame(info)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/positive_region_annotation.txt"
write.table(info,myoutf,sep="\t",quote=F)


[1]HepG2 TF signal
rm(list=ls())

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_integration_prediction.Rda")
res1 = data

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_TF_integration_Core.Rda")

com = intersect(colnames(res1),colnames(res2))
res2 = res2[,com]

tmpinf =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/positive_region_annotation.txt"
info = read.table(tmpinf,sep="\t",quote=NULL)
#info$sta = ceiling((info$start+1)/bin.size)
row.names(info) = paste0(info$seqnames,"_",info$start,"_",info$end)
data = res2
com = intersect(row.names(info),row.names(data))
data= data[com,]
info = info[com,]

tag = data$tag == "1"
data = data[tag,]
data = data[,2:ncol(data)]

#delete the low quality regions
avg = apply(data,1,mean)
tag = avg > 1
data = data[tag,]

variation = apply(data,1, var)
tag = variation > 0.5
data = data[tag,]

#begin to re-overlap
com = intersect(row.names(info),row.names(data))
info = info[com,]
data = data[com,]

#get the promoter
tag = grep("Promoter",info$annotation)
sam1 = row.names(info)[tag]
promoter = data[sam1,]

tmp = info[-tag,]
tag = grep("Exon",tmp$annotation)
tmp = tmp[-tag,]
sam2 = row.names(tmp)
distal = data[sam2,]

info = cbind(Name = c(row.names(promoter),row.names(distal)), Annotation = c(rep("Promoter",nrow(promoter)), rep("Distal",nrow(distal))))
info = as.data.frame(info)

data = rbind(promoter,distal)
raw.data = data

data = raw.data
tag = c(rep(1,nrow(promoter)),rep(0,nrow(distal)))
data = cbind(tag,data)
data$tag = factor(data$tag, levels =c(0,1))

#separate the positive and negative
tag1 = data$tag == "1"
tag2 = data$tag == "0"

dat1 = data[tag1,]
dat2 = data[tag2,]

tag1 = sample(nrow(dat1), 2000)
tag2 = sample(nrow(dat2), 2000)

data = rbind(dat1[tag1,],dat2[tag2,])
source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

TF_fit = ROC(data) 

#change to histone names
#change the name 
com = intersect(colnames(data), row.names(final_info))
data = data[,com]
final_info = final_info[com,]
colnames(data) = final_info$Experiment.target
tag = c(rep(1,2000),rep(0,2000))
data = cbind(tag,data)
data$tag = factor(data$tag, levels=c(0,1))

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_prediction_model.Rda")

res = NULL
tmp = predict(model, data[,-1], type="prob")
tmp = data.frame(data[,1], tmp)
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
	
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_prediction_promoter_distal.Rda"
save.image(file=myoutf)

#examine if features have issues
raw.data = data
data = data[,2:ncol(data)]

#principal component analysis
res = prcomp(data,scale=T)
res = res$x
res = res[,c(1,2)]


res = as.data.frame(res)
row.names(info) = info$Name

com = intersect(row.names(info),row.names(res))
info = info[com,]
res = res[com,]
res$Annotation = info$Annotation
colnames(res)[1:2] = c("x","y")

library(RColorBrewer)
c1 <- brewer.pal(10,"Set1")[2]
#c2 <- brewer.pal(10,"Set1")[2]
c2 <- "#999999"

px = pretty(res$x)
py = pretty(res$y)
commonTheme = list(labs(color="Density",fill="Density",
                        x="PC1",
                        y="PC2"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/HepG2/HepG2_TF_Proximal_vs_distal_association.pdf"
pdf(myoutf1, width= 3, height= 3)
p <- ggplot(res, aes(x=x,y=y,color=Annotation)) +   theme(legend.position="none")  +  geom_point(size=0.25)+scale_color_manual(values=c(c1,c2))
p1 <- p +  ggtitle("HepG2 \n (Proximal vs Distal)")
p2 <- p1 + theme(panel.grid.minor = element_line(colour="green", size=0.5),panel.grid.major = element_line(colour="green", size=0.5)) + theme_bw()
p3 <- p2 + guides(alpha="none",fill=F) + commonTheme + theme(aspect.ratio=1,panel.border = element_rect(colour = "black",size=1),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=1),axis.text=element_text(size= 7.25),legend.position="none") 
p3 <- p3 +  scale_x_continuous(breaks=px, limits=range(px)) + scale_y_continuous(breaks=py, limits=range(py))
#p3 <- p3 + annotate("text",x=px[length(px)-1],y=py[length(py)-1],label=paste0("rho=",coef),size=1.5,)
print(p3)
dev.off()

[2]Examine the HepG2 histone 
rm(list=ls())

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_integration_prediction.Rda")
res2 = data

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_Histone_integration_Core.Rda")

com = intersect(colnames(res1),colnames(res2))
res2 = res1[,com]

tmpinf =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/positive_region_annotation.txt"
info = read.table(tmpinf,sep="\t",quote=NULL)
#info$sta = ceiling((info$start+1)/bin.size)
row.names(info) = paste0(info$seqnames,"_",info$start,"_",info$end)
data = res2
com = intersect(row.names(info),row.names(data))
data= data[com,]
info = info[com,]

tag = data$tag == "1"
data = data[tag,]
data = data[,2:ncol(data)]

#delete the low quality regions
avg = apply(data,1,mean)
tag = avg > 1
data = data[tag,]

variation = apply(data,1, var)
tag = variation > 0.5
data = data[tag,]

#begin to re-overlap
com = intersect(row.names(info),row.names(data))
info = info[com,]
data = data[com,]

#get the promoter
tag = grep("Promoter",info$annotation)
sam1 = row.names(info)[tag]
promoter = data[sam1,]

tmp = info[-tag,]
tag = grep("Exon",tmp$annotation)
tmp = tmp[-tag,]
sam2 = row.names(tmp)
distal = data[sam2,]

info = cbind(Name = c(row.names(promoter),row.names(distal)), Annotation = c(rep("Promoter",nrow(promoter)), rep("Distal",nrow(distal))))
info = as.data.frame(info)

data = rbind(promoter,distal)
raw.data = data

data = raw.data
tag = c(rep(1,nrow(promoter)),rep(0,nrow(distal)))
data = cbind(tag,data)
data$tag = factor(data$tag, levels =c(0,1))

#separate the positive and negative
tag1 = data$tag == "1"
tag2 = data$tag == "0"

dat1 = data[tag1,]
dat2 = data[tag2,]

tag1 = sample(nrow(dat1), 2000)
tag2 = sample(nrow(dat2), 2000)

data = rbind(dat1[tag1,],dat2[tag2,])
source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")

HM_fit = ROC(data) 

#change to histone names
#change the name 
com = intersect(colnames(data), row.names(final_info))
data = data[,com]
final_info = final_info[com,]
colnames(data) = final_info$Experiment.target
tag = c(rep(1,2000),rep(0,2000))
data = cbind(tag,data)
data$tag = factor(data$tag, levels=c(0,1))

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_prediction_model.Rda")

res = NULL
tmp = predict(model, data[,-1], type="prob")
tmp = data.frame(data[,1], tmp)
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
	
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_prediction_promoter_distal.Rda"
save.image(file=myoutf)

#principal component analysis
data = data[,2:ncol(data)]
res = prcomp(data,scale = T)
res = res$x
res = res[,c(1,2)]

res = data.umap$layout
res = as.data.frame(res)
row.names(info) = info$Name

com = intersect(row.names(info),row.names(res))
info = info[com,]
res = res[com,]
res$Annotation = info$Annotation
colnames(res)[1:2] = c("x","y")

library(RColorBrewer)
c1 <- brewer.pal(10,"Set1")[1]
#c2 <- brewer.pal(10,"Set1")[2]
c2 <- "#999999"

px = pretty(res$x)
py = pretty(res$y)
commonTheme = list(labs(color="Density",fill="Density",
                        x="PC1",
                        y="PC2"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/HepG2/HepG2_Histone_Proximal_vs_distal_association.pdf"
pdf(myoutf1, width= 3, height= 3)
p <- ggplot(res, aes(x=x,y=y,color=Annotation)) +   theme(legend.position="none")  +  geom_point(size=0.25)+scale_color_manual(values=c(c1,c2))
p1 <- p +  ggtitle("HepG2 \n (Proximal vs Distal)")
p2 <- p1 + theme(panel.grid.minor = element_line(colour="green", size=0.5),panel.grid.major = element_line(colour="green", size=0.5)) + theme_bw()
p3 <- p2 + guides(alpha="none",fill=F) + commonTheme + theme(aspect.ratio=1,panel.border = element_rect(colour = "black",size=1),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=1),axis.text=element_text(size= 7.25),legend.position="none") 
p3 <- p3 +  scale_x_continuous(breaks=px, limits=range(px)) + scale_y_continuous(breaks=py, limits=range(py))
#p3 <- p3 + annotate("text",x=px[length(px)-1],y=py[length(py)-1],label=paste0("rho=",coef),size=1.5,)
print(p3)
dev.off()

[4]Begin to plot for the self model plot
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_prediction_promoter_distal.Rda"

load(myinf1)

xx1 = fit[[1]]
yy1 = fit[[2]]

TF_xx1 = 1-xx1
TF_yy1 = 1-yy1

auc1 = paste0("TF(AUC=", round(mean(1-fit[[3]]),2),")")

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_prediction_promoter_distal.Rda"
load(myinf2)
xx2 = fit[[1]]
yy2 = fit[[2]]

histone_xx2 = 1 - xx2
histone_yy2 = 1 - yy2

auc2 = paste0("Histone(AUC=", round(mean(1-fit[[3]]),2),")")
label = c(auc1,auc2)

color <- brewer.pal(10,"Set1")[1:2]
myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/HepG2/ROC_Curve_TF_histone_self_prediction.pdf"
pdf(myoutf, width= 5, height= 5)
par(pty="s")
plot(TF_xx1,TF_yy1,main="", xlab="False Positive Rate (1-Specificity)",ylab="True Positive Rate (Sensitivity)",
xlim=c(0,1),ylim=c(0,1),cex=0,cex.main=1.5,cex.lab=1,font.main=2) #xaxs="i",yaxs="i"
par(new=TRUE)
lines(TF_xx1,TF_yy1,lwd=2,col=color[1],lty=1)
par(new=TRUE)
lines(histone_xx2,histone_yy2,lwd=2,col=color[2],lty=2)
abline(0,1,lty=3)
legend("bottomright",label,lty=c(1,2),lwd=rep(2.5,length(label)),col=color[1:2],cex=0.75, box.lty=0)
dev.off()

[5]Change to barplot 
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_prediction_promoter_distal.Rda"
load(myinf1)
TF_auc = mean(1-TF_fit[[3]])
TF_auc_variation = var(1-TF_fit[[3]])

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_prediction_promoter_distal.Rda"
load(myinf2)
histone_auc = mean(1-HM_fit[[3]])
histone_auc_variation = var(1-HM_fit[[3]])

data = rbind(cbind(TF_auc,TF_auc_variation),cbind(histone_auc,histone_auc_variation))

data = as.data.frame(data)
data$target = c("TF","Histone")
colnames(data)[1:2] = c("AUC","Sd")
data$AUC = round(data$AUC,2)

#c1 <- brewer.pal(10,"RdBu")[3]

c1 <- brewer.pal(10,"Set1")[1]
c2 <- brewer.pal(10,"Set1")[2]
color <- c(c1,c2)
data$target = factor(data$target,levels=c("TF","Histone"))

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/HepG2/TF_histone_self_proximal_distal_prediction_barplot.pdf"
pdf(myoutf1, width= 3, height= 3)
p <- position_dodge(0.1)
p1 <- ggplot(data=data, aes(x=target, y=AUC, fill=color)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5) + geom_errorbar(aes(ymin=AUC-Sd, ymax=AUC+Sd), width=.2,position=position_dodge(.9)) + scale_fill_manual(values=c(c2,c1))+guides(fill=FALSE)+scale_y_continuous(expand = c(0,0),limits = c(0, 0.90)) + geom_text(aes(label=AUC), vjust=1.6, color="white", size=3.5)
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(angle=45,size= 7.25),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank()) + ylab("AUC score") + geom_hline(yintercept=0.5, linetype="dashed", color = "red", size=0.5)
p10 <- p9 +ggtitle(" ") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
p10 <- p10 + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
p10
dev.off()


[6]Begin to plot for the ATAC-seq model plot
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_prediction_promoter_distal.Rda"
load(myinf1)

TF_xx1 = 1- xx
TF_yy1 = 1- yy
auc1 = paste0("TF(AUC=", round(1-sum(tmp1*tmp2),2),")")

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_prediction_promoter_distal.Rda"
load(myinf2)
histone_xx2 = 1-xx
histone_yy2 = 1-yy
auc2 = paste0("Histone(AUC=", round(1-sum(tmp1*tmp2),2),")")
label = c(auc1,auc2)

color <- brewer.pal(10,"Set1")[1:2]

myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/HepG2/ROC_Curve_TF_histone_model_prediction.pdf"
pdf(myoutf, width= 5, height= 5)
par(pty="s")
plot(TF_xx1,TF_yy1,main="", xlab="False Positive Rate (1-Specificity)",ylab="True Positive Rate (Sensitivity)",
xlim=c(0,1),ylim=c(0,1),cex=0,cex.main=1.5,cex.lab=1,font.main=2) #xaxs="i",yaxs="i"
par(new=TRUE)
lines(TF_xx1,TF_yy1,lwd=2,col=color[1],lty=1)
par(new=TRUE)
lines(histone_xx2,histone_yy2,lwd=2,col=color[2],lty=2)
abline(0,1,lty=3)
legend("bottomright",label,lty=c(1,2),lwd=rep(2.5,length(label)),col=color[1:2],cex=0.75, box.lty=0)
dev.off()

[7]Begin to plot for other model in barplot
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_prediction_promoter_distal.Rda"
load(myinf1)
TF_auc = 1- sum(tmp1*tmp2)

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_prediction_promoter_distal.Rda"
load(myinf2)
histone_auc = 1- sum(tmp1*tmp2)

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

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/HepG2/TF_histone_other_Model_proximal_distal_prediction_barplot.pdf"
pdf(myoutf1, width= 3, height= 3)
p <- position_dodge(0.1)
p1 <- ggplot(data=data, aes(x=target, y=AUC, fill=color)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5) + scale_fill_manual(values=c(c2,c1))+guides(fill=FALSE)+scale_y_continuous(expand = c(0,0),limits = c(0, 0.90)) + geom_text(aes(label=AUC), vjust=1.6, color="white", size=3.5)
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(angle=45,size= 7.25),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank()) + ylab("AUC score") + geom_hline(yintercept=0.5, linetype="dashed", color = "red", size=0.5)
p10 <- p9 +ggtitle(" ") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
p10 <- p10 + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
p10
dev.off()


