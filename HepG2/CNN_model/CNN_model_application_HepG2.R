[1]Model training
rm(list=ls())
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/CNN/HepG2_local_training_CNN.Rda")
#install.packages("tensorflow")
library(tensorflow)
#install_tensorflow()

library(abind)
library(keras)
library("Biostrings")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model <- keras_model_sequential() 
model %>% 
  layer_conv_1d(filter = 320, kernel_size = 8, activation = 'relu', input_shape = c(100, 4)) %>% 
  layer_conv_1d(filter = 320, kernel_size = 8, activation = 'relu', input_shape = c(100, 4)) %>% 
  layer_dropout(rate = 0.2) %>%
  layer_max_pooling_1d(pool_size = 2) %>%
  layer_conv_1d(filter = 480, kernel_size = 8, activation = 'relu', input_shape = c(100, 4)) %>% 
  layer_conv_1d(filter = 480, kernel_size = 8, activation = 'relu', input_shape = c(100, 4)) %>% 
  layer_dropout(rate = 0.2) %>%
  layer_max_pooling_1d(pool_size = 2) %>%
  layer_conv_1d(filter = 640, kernel_size = 8, activation = 'relu', input_shape = c(100, 4)) %>%  
  layer_flatten() %>%
  layer_dropout(rate = 0.2) %>% 
  layer_dense(2000,activation = 'relu') %>% 
  layer_dense(120,activation = 'relu') %>% 
#  layer_dense(units = 2, activation = 'softmax')		## change units=1 and activation="sigmoid" for two-class classification
  layer_dense(units = 1, activation = 'sigmoid')

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#max pooling
#kernel size
#epoach
summary(model)

model %>% compile(
#  loss = 'categorical_crossentropy',							## loss = "binary_crossentropy" for for two-class classification
  loss = 'binary_crossentropy',							
  optimizer = optimizer_rmsprop(lr = 0.0001, decay = 1e-6),
  metrics = c('accuracy')
)


history <- model %>% fit(
  tr.x, tr.y, 
  epochs = 50, batch_size = 128, 
  validation_split = 0.2
)

myoutf ="/lorax/chenglab/yanding/ATAC_seq_integration/Figures/CNN/CNN_HepG2_training.pdf"
pdf(width=5,height=5,myoutf)
plot(history)
dev.off()

cla= model %>% predict_classes(te.x)
sum(cla==te.y)
sum(cla!=te.y)


prob = model %>% predict_proba(te.x)

xx = prob
cor1 = xx[te.y==1]
cor2 = xx[te.y==0]
length(cor1)
length(cor2)


thr = sort(xx)
yy =  xx =  rep(0, length(thr))
for(i in 1:length(thr))
{
	aa = sum(cor1>=thr[i])
	bb = sum(cor1<thr[i])
	cc = sum(cor2>=thr[i])
	dd = sum(cor2<thr[i])
	yy[i] = aa/(aa+bb)
	xx[i] = cc/(cc+dd)
}
xx = c(1, xx, 0)
yy = c(1, yy, 0)
tmp1 = tmp2 = rep(0,length(xx)-1)
for(i in 1:(length(xx)-1))
{
  tmp1[i] = xx[i]-xx[i+1]
  tmp2[i] = (yy[i+1]+yy[i])/2	
}
myauc = sum(tmp1*tmp2)
myauc

#begin to plot
label = c("AUC=0.63")


color <- brewer.pal(10,"Set1")[3]
myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_2/ROC_Curve_HepG2_CNN_training.pdf"
pdf(myoutf, width= 5, height= 5)
par(pty="s")
plot(xx,yy,main="", xlab="False Positive Rate (1-Specificity)",ylab="True Positive Rate (Sensitivity)",
xlim=c(0,1),ylim=c(0,1),cex=0,cex.main=1.5,cex.lab=1,font.main=2) #xaxs="i",yaxs="i"
par(new=TRUE)
lines(xx,yy,lwd=2,col=color[1],lty=1)
abline(0,1,lty=3)
legend("bottomright",label,lty=c(1),lwd=rep(2.5,length(label)),col=color[1],cex=0.75, box.lty=0)
dev.off()

#begin to apply the validation cohort
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/CNN/HepG2_local_validation_CNN.Rda")

prob = model %>% predict_proba(te.x)

xx = prob
cor1 = xx[te.y==1]
cor2 = xx[te.y==0]
length(cor1)
length(cor2)


thr = sort(xx)
yy =  xx =  rep(0, length(thr))
for(i in 1:length(thr))
{
	aa = sum(cor1>=thr[i])
	bb = sum(cor1<thr[i])
	cc = sum(cor2>=thr[i])
	dd = sum(cor2<thr[i])
	yy[i] = aa/(aa+bb)
	xx[i] = cc/(cc+dd)
}
xx = c(1, xx, 0)
yy = c(1, yy, 0)
tmp1 = tmp2 = rep(0,length(xx)-1)
for(i in 1:(length(xx)-1))
{
  tmp1[i] = xx[i]-xx[i+1]
  tmp2[i] = (yy[i+1]+yy[i])/2	
}
myauc = sum(tmp1*tmp2)
myauc

#begin to plot
label = c("AUC=0.63")


color <- brewer.pal(10,"Set1")[3]
myoutf <- "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_2/ROC_Curve_HepG2_CNN_validation.pdf"
pdf(myoutf, width= 5, height= 5)
par(pty="s")
plot(xx,yy,main="", xlab="False Positive Rate (1-Specificity)",ylab="True Positive Rate (Sensitivity)",
xlim=c(0,1),ylim=c(0,1),cex=0,cex.main=1.5,cex.lab=1,font.main=2) #xaxs="i",yaxs="i"
par(new=TRUE)
lines(xx,yy,lwd=2,col=color[1],lty=1)
abline(0,1,lty=3)
legend("bottomright",label,lty=c(1),lwd=rep(2.5,length(label)),col=color[1],cex=0.75, box.lty=0)
dev.off()

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/CNN/HepG2_local_CNN_model.Rda"
save.image(myoutf)