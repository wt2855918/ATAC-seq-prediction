[1]GM12878 co-localization
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda"
load(myinf1)
com = intersect(colnames(data),row.names(final_info))
data = data[,com]
final_info = final_info[com,]
colnames(data) = final_info$Experiment.target
TF_data = data

#alternative TF data
#myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_TF_data.Rda"
#load(myinf1)
#TF_data = GM12878_data

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_integration_prediction.Rda"
load(myinf2)
Histone_data = data
Histone_data = Histone_data[,2:ncol(Histone_data)]
colnames(Histone_data) = final_info$Experiment.target


com = intersect(row.names(TF_data),row.names(Histone_data))
TF_data = TF_data[com,]
Histone_data = Histone_data[com,]

#########Pheatmap####################
res = matrix(0,ncol(Histone_data),ncol(TF_data))
row.names(res) = colnames(Histone_data)
colnames(res) = colnames(TF_data)
res = as.data.frame(res)

for(i in 1 : ncol(TF_data))
{
	cat("\r",i)
	
	xx1 = TF_data[,i]
	fit = cor(xx1,Histone_data)

	res[,i] = as.numeric(fit)
}

res = as.matrix(res)

library(pheatmap)


data = as.matrix(data)
diag(data) = NA

library(RColorBrewer)
c1 <- brewer.pal(6,"YlOrRd")
c1 <- brewer.pal(8,"YlGnBu")
c1 <- c(brewer.pal(6,"YlGnBu"),brewer.pal(8,"PuBu")[7:8])

c1 <- brewer.pal(6,"YlGnBu")


c1 <- brewer.pal(8,"YlGnBu")[1:6]
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_6/TF_vs_Histon_correlation_GM12878.pdf"
pdf(myoutf,width=5.75,height=5)
pheatmap(res, cluster_rows = T, cluster_cols = T, show_rownames = F, show_colnames = F)
dev.off()

################Scatterplot###################
xx1 = apply(Histone_data,1,mean)
xx2 = apply(TF_data,1,mean)

xx1 = log2(xx1+1)
xx2 = log2(xx2+1)

df = cbind(xx1,xx2)
colnames(df) = c("x","y")
df = as.data.frame(df)

px <- pretty(df$x)
py <- pretty(df$y)

#get the coef
fit <- cor(df$x, df$y, method="s",use="complete.obs")
coef <- round(fit,2)

library(RColorBrewer)
c1 <- brewer.pal(10,"RdBu")[3]

#c1 <- "#F57E20"
commonTheme = list(labs(color="Density",fill="Density",
                        x="Histone",
                        y="TF"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))
c1 <- brewer.pal(10,"RdBu")[3]

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_6/GM12878_histone_TF_scatterplot.pdf"
pdf(myoutf1, width= 1.55, height= 1.55)
p <- ggplot(df, aes(x=x,y=y,color=c1)) +   theme(legend.position="none")  +  geom_point(size=0.05)+scale_color_manual(values=c(c1))
p1 <- p + geom_smooth(method=lm,linetype=1,colour="grey70",se=F,size=1) + ggtitle("GM12878 \n (TF vs Histone)")
p2 <- p1 + theme(panel.grid.minor = element_line(colour="green", size=0.5),panel.grid.major = element_line(colour="green", size=0.5)) + theme_bw()
p3 <- p2 + guides(alpha="none",fill=F) + commonTheme + theme(aspect.ratio=1,panel.border = element_rect(colour = "black",size=1),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=1),axis.text=element_text(size= 7.25),legend.position="none") 
p3 <- p3 +  scale_x_continuous(breaks=px, limits=range(px)) + scale_y_continuous(breaks=py, limits=range(py))
p3 <- p3 + annotate("text",x=px[length(px)-1],y=py[length(py)-1],label=paste0("rho=",coef),size=1.5,)
print(p3)
dev.off()


[2]GM12878 Histone + TF prediction
rm(list=ls())
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/TF_integration_prediction.Rda")
TF_AUC = 1 - fit[[3]]

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_integration_prediction.Rda")
Histone_AUC = 1 - fit[[3]]

myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Final_ID/Histone_TF_integration_prediction.Rda"
load(myinf3)
TF_HM_AUC = 1 - fit[[3]]

#relative importance

data = matrix(0,3,3)
row.names(data) = c("Histone","TF","Histone+TF")
colnames(data) = c("Target","Value","Sd")
data = as.data.frame(data)

row_1 = c("Histone",mean(Histone_AUC),sd(Histone_AUC))
row_2 = c("TF",mean(TF_AUC),sd(TF_AUC))
row_3 = c("Target",mean(TF_HM_AUC),sd(TF_HM_AUC))

data[1,] = row_1
data[2,] = row_2
data[3,] = row_3
data$Value = as.numeric(data$Value)
data$Sd = as.numeric(data$Sd)
data$Target = c("Histone","TF","Histone+TF")
data$Target = factor(data$Target, levels = c("Histone","TF","Histone+TF"))
data$Value = round(data$Value,2)

library(RColorBrewer)
c1 <- brewer.pal(10,"RdBu")[3]

myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_6/Histone_prediction_combined_GM12878.pdf"
pdf(myoutf, width= 3, height= 3)
p <- position_dodge(0.1)
p1 <- ggplot(data=data, aes(x=Target, y=Value, fill=c1)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5) + scale_fill_manual(values=c(c1,c1,c1))+guides(fill=FALSE)+scale_y_continuous(expand = c(0,0),limits = c(0, 0.865)) + geom_text(aes(label=Value), vjust=1.6, color="white", size=3.5)
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(angle=45,size= 7.25),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank()) + ylab("AUC score") + geom_hline(yintercept=0.5, linetype="dashed", color = "red", size=0.5)
p10 <- p9 +ggtitle(" ") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
p10 <- p10 + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
p10
dev.off()

[3]HepG2 scatterplot co-localization
rm(list=ls())
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_integration_prediction.Rda"
load(myinf1)
com = intersect(colnames(data),row.names(final_info))
data = data[,com]
final_info = final_info[com,]
colnames(data) = final_info$Experiment.target
TF_data = data

#alternative TF data
#myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Integration/GM12878_shared_TF_data.Rda"
#load(myinf1)
#TF_data = GM12878_data

myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_integration_prediction.Rda"
load(myinf2)
Histone_data = data
Histone_data = Histone_data[,2:ncol(Histone_data)]
colnames(Histone_data) = final_info$Experiment.target


com = intersect(row.names(TF_data),row.names(Histone_data))
TF_data = TF_data[com,]
Histone_data = Histone_data[com,]

#########Pheatmap####################
res = matrix(0,ncol(Histone_data),ncol(TF_data))
row.names(res) = colnames(Histone_data)
colnames(res) = colnames(TF_data)
res = as.data.frame(res)

for(i in 1 : ncol(TF_data))
{
	cat("\r",i)
	
	xx1 = TF_data[,i]
	fit = cor(xx1,Histone_data)

	res[,i] = as.numeric(fit)
}

res = as.matrix(res)

library(pheatmap)


data = as.matrix(data)
diag(data) = NA

library(RColorBrewer)
c1 <- brewer.pal(6,"YlOrRd")
c1 <- brewer.pal(8,"YlGnBu")
c1 <- c(brewer.pal(6,"YlGnBu"),brewer.pal(8,"PuBu")[7:8])

c1 <- brewer.pal(6,"YlGnBu")


c1 <- brewer.pal(8,"YlGnBu")[1:6]
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_6/TF_vs_Histon_correlation_HepG2.pdf"
pdf(myoutf,width=5.75,height=5)
pheatmap(res, cluster_rows = T, cluster_cols = T, show_rownames = F, show_colnames = F)
dev.off()

################Scatterplot###################
xx1 = apply(Histone_data,1,mean)
xx2 = apply(TF_data,1,mean)

xx1 = log2(xx1+1)
xx2 = log2(xx2+1)

df = cbind(xx1,xx2)
colnames(df) = c("x","y")
df = as.data.frame(df)

px <- pretty(df$x)
py <- pretty(df$y)

#get the coef
fit <- cor(df$x, df$y, method="s",use="complete.obs")
coef <- round(fit,2)

library(RColorBrewer)
c1 <- brewer.pal(10,"RdBu")[8]

#c1 <- "#F57E20"
commonTheme = list(labs(color="Density",fill="Density",
                        x="Histone",
                        y="TF"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))
c1 <- brewer.pal(10,"RdBu")[8]

myoutf1 <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_6/HepG2_histone_TF_scatterplot.pdf"
pdf(myoutf1, width= 1.55, height= 1.55)
p <- ggplot(df, aes(x=x,y=y,color=c1)) +   theme(legend.position="none")  +  geom_point(size=0.05)+scale_color_manual(values=c(c1))
p1 <- p + geom_smooth(method=lm,linetype=1,colour="grey70",se=F,size=1) + ggtitle("HepG2 \n (TF vs Histone)")
p2 <- p1 + theme(panel.grid.minor = element_line(colour="green", size=0.5),panel.grid.major = element_line(colour="green", size=0.5)) + theme_bw()
p3 <- p2 + guides(alpha="none",fill=F) + commonTheme + theme(aspect.ratio=1,panel.border = element_rect(colour = "black",size=1),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=1),axis.text=element_text(size= 7.25),legend.position="none") 
p3 <- p3 +  scale_x_continuous(breaks=px, limits=range(px)) + scale_y_continuous(breaks=py, limits=range(py))
p3 <- p3 + annotate("text",x=px[length(px)-1],y=py[length(py)-1],label=paste0("rho=",coef),size=1.5,)
print(p3)
dev.off()

[4]Barplot showing the same AUC HepG2
rm(list=ls())
rm(list=ls())
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/TF_integration_prediction.Rda")
TF_AUC = 1 - fit[[3]]

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_integration_prediction.Rda")
Histone_AUC = 1 - fit[[3]]

myinf3 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Final_ID/Histone_TF_integration_prediction.Rda"
load(myinf3)
TF_HM_AUC = 1 - fit[[3]]

#relative importance

data = matrix(0,3,3)
row.names(data) = c("Histone","TF","Histone+TF")
colnames(data) = c("Target","Value","Sd")
data = as.data.frame(data)

row_1 = c("Histone",mean(Histone_AUC),sd(Histone_AUC))
row_2 = c("TF",mean(TF_AUC),sd(TF_AUC))
row_3 = c("Target",mean(TF_HM_AUC),sd(TF_HM_AUC))

data[1,] = row_1
data[2,] = row_2
data[3,] = row_3
data$Value = as.numeric(data$Value)
data$Sd = as.numeric(data$Sd)
data$Target = c("Histone","TF","Histone+TF")
data$Target = factor(data$Target, levels = c("Histone","TF","Histone+TF"))
data$Value = round(data$Value,2)

library(RColorBrewer)
c1 <- brewer.pal(10,"RdBu")[8]

myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_6/Histone_prediction_combined_HepG2.pdf"
pdf(myoutf, width= 3, height= 3)
p <- position_dodge(0.1)
p1 <- ggplot(data=data, aes(x=Target, y=Value, fill=c1)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5) + scale_fill_manual(values=c(c1,c1,c1))+guides(fill=FALSE)+scale_y_continuous(expand = c(0,0),limits = c(0, 0.865)) + geom_text(aes(label=Value), vjust=1.6, color="white", size=3.5)
p9 <- p1 + theme(axis.text.y=element_text(size= 7.25),axis.text.x=element_text(angle=45,size= 7.25),axis.title.y=element_text(size= 7.25),axis.title.x=element_blank()) + ylab("AUC score") + geom_hline(yintercept=0.5, linetype="dashed", color = "red", size=0.5)
p10 <- p9 +ggtitle(" ") + theme(panel.border = element_blank(),axis.line = element_line(colour="black"),axis.title=element_text(size=7.25),plot.title=element_text(size=7.25,hjust = 0.5,vjust=0.5))
p10 <- p10 + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
p10
dev.off()
