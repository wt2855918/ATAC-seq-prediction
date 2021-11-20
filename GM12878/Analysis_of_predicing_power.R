[1]Histone modification
rm(list=ls())
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_integration_prediction.Rda"
load(myinf)

tmpinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/Histone_file_annotation_bigwig.txt"
#info = read.table(tmpinf,sep="\t",quote=NULL)

#output the AUC with standard bar

#AUC curve cell line
label <- paste0("AUC=",round(mean(1-fit[[3]]),2))
color <- brewer.pal(10,"Set1")[1]
myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/Histone_prediction.pdf"
pdf(myoutf, width= 5, height= 5)
par(pty="s")
plot(1-fit[[1]],1-fit[[2]],main="Histone", xlab="False Positive Rate (1-Specificity)",ylab="True Positive Rate (Sensitivity)",
xlim=c(0,1),ylim=c(0,1),cex=0,cex.main=1.5,cex.lab=1,font.main=2) #xaxs="i",yaxs="i"
par(new=TRUE)
lines(1-fit[[1]],1-fit[[2]],lwd=2,col=color[1],lty=1)
abline(0,1,lty=3)
legend("bottomright",label,lty=rep(1,1),lwd=rep(2.5,1),col=color[1],cex=0.75, box.lty=0)
dev.off()

#relative importance
res = fit[[6]]
info = info[row.names(res),]
row.names(res) = info$target

data = matrix(0,nrow(res),3)
row.names(data) = row.names(res)
colnames(data) = c("Target","Value","Sd")
data = as.data.frame(data)

for(i in 1 : nrow(res))
{
	cat("\r",i)
	
	xx = res[i,]
	data[i,"Value"] = mean(xx)
	data[i,"Sd"] = sd(xx)
	data[i,"Target"] = row.names(res)[i]
	
}
data = data[order(data$Value,decreasing=T),]
data$ID = row.names(data)
data$ID = factor(data$ID, levels = data$ID)

myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/Histone_prediction_importance.pdf"
pdf(myoutf, width= 7.5, height= 4.5)

p <- ggplot(data, aes(x=ID, y=Value, fill=Target)) + 
   geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Value-Sd, ymax=Value+Sd), width=.2,
                 position=position_dodge(.9))
p <- p + theme_minimal() + theme(axis.text.x=element_text(angle=90, hjust=1))
p

dev.off()

[2]TF modification
rm(list=ls())
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda"
load(myinf)

tmpinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/TF_file_annotation_bigwig.txt"
#info = read.table(tmpinf,sep="\t",quote=NULL)

#output the AUC with standard bar

#AUC curve cell line
label <- paste0("AUC=",round(mean(1-fit[[3]]),2))
color <- brewer.pal(10,"Set1")[2]
myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/TF_prediction.pdf"
pdf(myoutf, width= 5, height= 5)
par(pty="s")
plot(1-fit[[1]],1-fit[[2]],main="TF", xlab="False Positive Rate (1-Specificity)",ylab="True Positive Rate (Sensitivity)",
xlim=c(0,1),ylim=c(0,1),cex=0,cex.main=1.5,cex.lab=1,font.main=2) #xaxs="i",yaxs="i"
par(new=TRUE)
lines(1-fit[[1]],1-fit[[2]],lwd=2,col=color[1],lty=1)
abline(0,1,lty=3)
legend("bottomright",label,lty=rep(1,1),lwd=rep(2.5,1),col=color[1],cex=0.75, box.lty=0)
dev.off()

#relative importance
res = fit[[6]]
info = info[row.names(res),]
row.names(res) = info$target

data = matrix(0,nrow(res),3)
row.names(data) = row.names(res)
colnames(data) = c("Target","Value","Sd")
data = as.data.frame(data)

for(i in 1 : nrow(res))
{
	cat("\r",i)
	
	xx = res[i,]
	data[i,"Value"] = mean(xx)
	data[i,"Sd"] = sd(xx)
	data[i,"Target"] = row.names(res)[i]
	
}
data = data[order(data$Value,decreasing=T),]
data$ID = row.names(data)
data = data[1:20,]

data$ID = factor(data$ID, levels = data$ID)

myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/TF_prediction_importance.pdf"
pdf(myoutf, width= 7.5, height= 4.5)

p <- ggplot(data, aes(x=ID, y=Value, fill=Target)) + 
   geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Value-Sd, ymax=Value+Sd), width=.2,
                 position=position_dodge(.9))
p <- p + theme_minimal() + theme(axis.text.x=element_text(angle=90, hjust=1))
p
dev.off()

[3]Motif prediction
rm(list=ls())
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/motif_integration_prediction.Rda"
load(myinf)

tmpinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/motif_file_annotation_bigwig.txt"

#output the AUC with standard bar

#AUC curve cell line
label <- paste0("AUC=",round(mean(1-fit[[3]]),2))
color <- brewer.pal(10,"Set1")[2]
myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/motif_prediction.pdf"
pdf(myoutf, width= 5, height= 5)
par(pty="s")
plot(1-fit[[1]],1-fit[[2]],main="motif", xlab="False Positive Rate (1-Specificity)",ylab="True Positive Rate (Sensitivity)",
xlim=c(0,1),ylim=c(0,1),cex=0,cex.main=1.5,cex.lab=1,font.main=2) #xaxs="i",yaxs="i"
par(new=TRUE)
lines(1-fit[[1]],1-fit[[2]],lwd=2,col=color[1],lty=1)
abline(0,1,lty=3)
legend("bottomright",label,lty=rep(1,1),lwd=rep(2.5,1),col=color[1],cex=0.75, box.lty=0)
dev.off()

#relative importance
res = fit[[6]]
row.names(res) = colnames(pos_count)

data = matrix(0,nrow(res),3)
row.names(data) = row.names(res)
colnames(data) = c("Target","Value","Sd")
data = as.data.frame(data)

for(i in 1 : nrow(res))
{
	cat("\r",i)
	
	xx = res[i,]
	data[i,"Value"] = mean(xx)
	data[i,"Sd"] = sd(xx)
	data[i,"Target"] = row.names(res)[i]
	
}
data = data[order(data$Value,decreasing=T),]
data$ID = row.names(data)
data = data[1:20,]
data$ID = factor(data$ID, levels = data$ID)

myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/motif_prediction_importance.pdf"
pdf(myoutf, width= 7.5, height= 4.5)

p <- ggplot(data, aes(x=ID, y=Value, fill=Target)) + 
   geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Value-Sd, ymax=Value+Sd), width=.2,
                 position=position_dodge(.9))
p <- p + theme_minimal() + theme(axis.text.x=element_text(angle=90, hjust=1))
p
dev.off()

[4]sequence prediction
rm(list=ls())
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Sequence_integration_prediction.Rda")

#AUC curve cell line
label <- paste0("AUC=",round(mean(1-fit[[3]]),2))
color <- brewer.pal(10,"Set1")[1]
myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/sequence_prediction.pdf"
pdf(myoutf, width= 5, height= 5)
par(pty="s")
plot(1-fit[[1]],1-fit[[2]],main="Sequence", xlab="False Positive Rate (1-Specificity)",ylab="True Positive Rate (Sensitivity)",
xlim=c(0,1),ylim=c(0,1),cex=0,cex.main=1.5,cex.lab=1,font.main=2) #xaxs="i",yaxs="i"
par(new=TRUE)
lines(1-fit[[1]],1-fit[[2]],lwd=2,col=color[1],lty=1)
abline(0,1,lty=3)
legend("bottomright",label,lty=rep(1,1),lwd=rep(2.5,1),col=color[1],cex=0.75, box.lty=0)
dev.off()

[5]Conservation score
rm(list=ls())
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/conScore_integration_prediction.Rda")

#AUC curve cell line
label <- paste0("AUC=",round(mean(1-fit[[3]]),2))
color <- brewer.pal(10,"Set1")[2]
myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/conScore_prediction.pdf"
pdf(myoutf, width= 5, height= 5)
par(pty="s")
plot(1-fit[[1]],1-fit[[2]],main="conScore", xlab="False Positive Rate (1-Specificity)",ylab="True Positive Rate (Sensitivity)",
xlim=c(0,1),ylim=c(0,1),cex=0,cex.main=1.5,cex.lab=1,font.main=2) #xaxs="i",yaxs="i"
par(new=TRUE)
lines(1-fit[[1]],1-fit[[2]],lwd=2,col=color[1],lty=1)
abline(0,1,lty=3)
legend("bottomright",label,lty=rep(1,1),lwd=rep(2.5,1),col=color[1],cex=0.75, box.lty=0)
dev.off()

[5]Integrative model
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_TF_conScore_integration_prediction.Rda"
load(myinf1)
model1 = 1-fit[[3]]

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_TF_conScore_sequence_integration_prediction.Rda"
load(myinf2)
model2 = 1-fit[[3]]

myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Without_seq_integration_prediction.Rda"
load(myinf3)
model3 = 1 -fit[[3]]

myinf4 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Final_integration_prediction.Rda"
load(myinf4)
model4 = 1 - fit[[3]]

xx1 = round(mean(model1),3)
sd1 = sd(model1)

xx2 = round(mean(model2),3)
sd2 = sd(model2)

xx3 = round(mean(model3),3)
sd3 = sd(model3)

xx4 = round(mean(model4),3)
sd4 = sd(model4)

data = matrix(0,4,2)
colnames(data) = c("AUC","Sd")
data = as.data.frame(data)

data$AUC = c(xx1,xx2,xx3,xx4)
data$Sd = c(sd1,sd2,sd3,sd4)

data$ID = paste0("Model_",row.names(data))
data$ID = factor(data$ID, levels = data$ID)

myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/integrative_prediction_AUC.pdf"
pdf(myoutf, width= 3, height= 3)

p <- ggplot(data, aes(x=ID, y=AUC)) + 
   geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=AUC-Sd, ymax=AUC+Sd), width=.2,
                 position=position_dodge(.9))
p <- p + theme_minimal() + theme(axis.text.x=element_text(angle=90, hjust=1)) + coord_cartesian(ylim=c(0.8, 0.88))
p
dev.off()

[6]Begin to annotate the peak 
rm(list=ls())

#induce the function 
#define the convert function
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


library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_positive.Rda")

positive$V2 = as.numeric(as.vector(positive$V2))
positive$V3 = as.numeric(as.vector(positive$V3))

peak = peakDF2GRanges(positive)
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
slotNames(peakAnno)
res = peakAnno@anno
res = as.data.frame(res)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_positive_annotation.txt"
write.table(res,myoutf,sep="\t",quote=F)


[7]Different signal difference

load ("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Final_integration_prediction_promoter_distal.Rda")
#AUC curve cell line
label <- paste0("AUC=",round(mean(1-fit[[3]]),2))
color <- brewer.pal(10,"Set1")[2]
myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/promoter_vs_non_promoter_prediction.pdf"
pdf(myoutf, width= 5, height= 5)
par(pty="s")
plot(1-fit[[1]],1-fit[[2]],main="Model4", xlab="False Positive Rate (1-Specificity)",ylab="True Positive Rate (Sensitivity)",
xlim=c(0,1),ylim=c(0,1),cex=0,cex.main=1.5,cex.lab=1,font.main=2) #xaxs="i",yaxs="i"
par(new=TRUE)
lines(1-fit[[1]],1-fit[[2]],lwd=2,col=color[1],lty=1)
abline(0,1,lty=3)
legend("bottomright",label,lty=rep(1,1),lwd=rep(2.5,1),col=color[1],cex=0.75, box.lty=0)
dev.off()

#relative importance
res = fit[[6]]

data = matrix(0,nrow(res),3)
row.names(data) = row.names(res)
colnames(data) = c("Target","Value","Sd")
data = as.data.frame(data)

for(i in 1 : nrow(res))
{
	cat("\r",i)
	
	xx = res[i,]
	data[i,"Value"] = mean(xx)
	data[i,"Sd"] = sd(xx)
	data[i,"Target"] = row.names(res)[i]
	
}
data = data[order(data$Value,decreasing=T),]
data$ID = row.names(data)
data = data[1:20,]
#data$ID = factor(data$ID, levels = data$ID)

tmpinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/Histone_file_annotation_bigwig.txt"
tmpinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/TF_file_annotation_bigWig.txt"

info1 = read.table(tmpinf1,sep="\t",quote=NULL)
info2 = read.table(tmpinf2,sep="\t",quote=NULL)

info1 = info1[,c("file_accession","target")]
info2 = info2[,c("file_accession","target")]
info = rbind(info1,info2)
row.names(info) = info$file_accession

com = intersect(row.names(data),row.names(info))
data[com,"ID"] = as.vector(info[com,"target"])
data$ID = factor(data$ID, levels = data$ID)
myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/promoter_vs_non_promoter_prediction_importance.pdf"
pdf(myoutf, width= 7.5, height= 4.5)

p <- ggplot(data, aes(x=ID, y=Value, fill=Target)) + 
   geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Value-Sd, ymax=Value+Sd), width=.2,
                 position=position_dodge(.9))
p <- p + theme_minimal() + theme(axis.text.x=element_text(angle=90, hjust=1))
p
dev.off()

[7]Begin to perform throw each histone out for power test
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

myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/Histone_file_annotation_bigwig.txt"
res = read.table(myinf,sep="\t",quote=NULL)

row.names(info) = info$sample_name
row.names(res) = res$file_accession

com = intersect(row.names(info),row.names(res))
info = info[com,]
res = res[com,]

info$Target = res$target
raw.info = info

load(tmpinf)
all_target = c("GM12878","0","All",mean(1-fit[[3]]),"All")

info = raw.info
info = apply(info,2, function(x) as.vector(x))
info = rbind(all_target, info)
info = as.data.frame(info)

row.names(info) = info$sample_name
info$Order_type = as.numeric(as.vector(info$Order_type))
info$AUC = as.numeric(as.vector(info$AUC))

info = info[order(info$Order_type),]
info$Name = paste0(info$sample_name,"_",info$Target)
info$Name = factor(info$Name, levels = info$Name)
info$AUC = round(info$AUC,2)
c1 <- brewer.pal(10,"RdBu")[3]

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/Histone_KO_each.pdf"
pdf(myoutf1, width= 10, height= 3)
p <- position_dodge(0.1)
p1 <- ggplot(data=info, aes(x=Name, y=AUC, fill=c1)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + scale_fill_manual(values=c(c1))+guides(fill=FALSE)+scale_y_continuous(expand = c(0,0),limits = c(0, 0.85)) + geom_text(aes(label=AUC), vjust=1.6, color="white", size=2)
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(angle=45,size= 7.25,hjust=1),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank()) + ylab("AUC score") + geom_hline(yintercept=0.5, linetype="dashed", color = "red", size=0.5)
p10 <- p9 +ggtitle(" ") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
p10 <- p10 + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
p10
dev.off()






