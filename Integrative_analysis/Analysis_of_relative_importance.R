[1]Barplot of GM12878 model
rm(list=ls())
myinf =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_Histone_integration_AUC.Rda"
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
colnames(info)[2] = "target"
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
p <- p + theme_minimal() + theme(axis.text.x=element_text(angle=90, hjust=1)) + ylab("Importance")
p

dev.off()

res = aggregate(Value ~ Target, data=data, mean)
res = res[order(res$Value,decreasing=T),]
     Target      Value
6   H3K4me2 17.7703345
1     H2AFZ 12.4953295
2   H3K27ac 11.1726283
5   H3K4me1  6.3625222
7   H3K4me3  6.3500751
8  H3K79me2  4.2491397
9    H3K9ac  3.1363406
10  H3K9me3  2.1146805
4  H3K36me3  1.4696000
3  H3K27me3  0.8431026
11 H4K20me1 -1.4738448

res$Target = factor(res$Target, levels = res$Target)

myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/GM12878/Histone_prediction_importance_mean.pdf"
pdf(myoutf, width= 7.5, height= 4.5)

p <- ggplot(res, aes(x=Target, y=Value, fill=Target)) + 
   geom_bar(stat="identity", position=position_dodge()) 
p <- p + theme_minimal() + theme(axis.text.x=element_text(angle=90, hjust=1)) + ylab("Importance")
p

dev.off()


[2]Barplot of HepG2 model
rm(list=ls())
myinf =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/HepG2_Histone_integration_AUC.Rda"
load(myinf)

tmpinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/Histone_file_annotation_bigwig.txt"
#info = read.table(tmpinf,sep="\t",quote=NULL)

#output the AUC with standard bar

#AUC curve cell line
label <- paste0("AUC=",round(mean(1-fit[[3]]),2))
color <- brewer.pal(10,"Set1")[1]
myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/HepG2/Histone_prediction.pdf"
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
colnames(info)[2] = "target"
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

myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/HepG2/Histone_prediction_importance.pdf"
pdf(myoutf, width= 7.5, height= 4.5)

p <- ggplot(data, aes(x=ID, y=Value, fill=Target)) + 
   geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Value-Sd, ymax=Value+Sd), width=.2,
                 position=position_dodge(.9))
p <- p + theme_minimal() + theme(axis.text.x=element_text(angle=90, hjust=1)) + ylab("Importance")
p

dev.off()

res = aggregate(Value ~ Target, data=data, mean)
res = res[order(res$Value,decreasing=T),]

     Target      Value
2   H3K27ac 11.9338703
6   H3K4me2 10.5465045
1     H2AFZ  9.9958460
7   H3K4me3  8.6409609
5   H3K4me1  7.6029602
9    H3K9ac  4.3421316
11 H4K20me1  3.9976983
8  H3K79me2  1.5877849
3  H3K27me3  1.1449371
10  H3K9me3  0.4958072
4  H3K36me3 -1.0303955

res$Target = factor(res$Target, levels = res$Target)

myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/HepG2/Histone_prediction_importance_mean.pdf"
pdf(myoutf, width= 7.5, height= 4.5)

p <- ggplot(res, aes(x=Target, y=Value, fill=Target)) + 
   geom_bar(stat="identity", position=position_dodge()) 
p <- p + theme_minimal() + theme(axis.text.x=element_text(angle=90, hjust=1)) + ylab("Importance")
p

dev.off()

[3]Barplot of mESC model
rm(list=ls())
myinf =  "/lorax/chenglab/yanding/ATAC_seq_integration/Data/mESC/mESC_Histone_integration_AUC.Rda"
load(myinf)

tmpinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/mESC/Histone_file_annotation_bigwig.txt"
#info = read.table(tmpinf,sep="\t",quote=NULL)

#output the AUC with standard bar

#AUC curve cell line
label <- paste0("AUC=",round(mean(1-fit[[3]]),2))
color <- brewer.pal(10,"Set1")[1]
myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/mESC/Histone_prediction.pdf"
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
colnames(info)[2] = "target"
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

myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/mESC/Histone_prediction_importance.pdf"
pdf(myoutf, width= 7.5, height= 4.5)

p <- ggplot(data, aes(x=ID, y=Value, fill=Target)) + 
   geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Value-Sd, ymax=Value+Sd), width=.2,
                 position=position_dodge(.9))
p <- p + theme_minimal() + theme(axis.text.x=element_text(angle=90, hjust=1)) + ylab("Importance")
p

dev.off()

res = aggregate(Value ~ Target, data=data, mean)
res = res[order(res$Value,decreasing=T),]

  Target    Value
4  H3K4me3 39.72393
2 H3K36me3 21.78856
5  H3K9me3 21.75092
1 H3K27me3 19.91544
3  H3K4me1 10.99064

res$Target = factor(res$Target, levels = res$Target)

myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/mESC/Histone_prediction_importance_mean.pdf"
pdf(myoutf, width= 7.5, height= 4.5)

p <- ggplot(res, aes(x=Target, y=Value, fill=Target)) + 
   geom_bar(stat="identity", position=position_dodge()) 
p <- p + theme_minimal() + theme(axis.text.x=element_text(angle=90, hjust=1)) + ylab("Importance")
p

dev.off()
