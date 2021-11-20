myinf1 ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1_training_and_validation/ENCFF172DEA_positive_training.fa"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1_training_and_validation/ENCFF172DEA_positive_validation.fa"
myinf3 ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1_training_and_validation/ENCFF172DEA_negative_training.fa"
myinf4 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1_training_and_validation/ENCFF172DEA_negative_validation.fa"

[1]build CNN

module load python/3.7-Anaconda
conda remove -name Velocity --all
source activate ATAC_seq_integration
conda install tensorflow

rm(list=ls())

#install.packages("tensorflow")
library(tensorflow)
#install_tensorflow()
library("Biostrings")

myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1_training_and_validation/ENCFF172DEA_positive_training.fa"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1_training_and_validation/ENCFF172DEA_random_negative_training.fa"

library(abind)
library(keras)
#---------------------------

data <- readDNAStringSet(myinf1)
nn = length(data)/2
data = as.data.frame(data)
tmp1 = row.names(data)
tmp2 = data$x
names(tmp2) = tmp1
data = tmp2
res = array(0, c(length(data), 100, 4))
for(k in 1:length(data))
{
	cat("\r", k)
	xx = unlist(strsplit(data[k], ""))
	se = which(xx=="A")
	res[k,se,1]=1
	se = which(xx=="T")
	res[k,se,2]=1
	se = which(xx=="C")
	res[k,se,3]=1
	se = which(xx=="G")
	res[k,se,4]=1	
}
dat.pos = res


data <- readDNAStringSet(myinf2)
nn = length(data)/2
data = as.data.frame(data)
tmp1 = row.names(data)
tmp2 = data$x
names(tmp2) = tmp1
data = tmp2
res = array(0, c(length(data), 100, 4))
for(k in 1:length(data))
{
	cat("\r", k)
	xx = unlist(strsplit(data[k], ""))
	se = which(xx=="A")
	res[k,se,1]=1
	se = which(xx=="T")
	res[k,se,2]=1
	se = which(xx=="C")
	res[k,se,3]=1
	se = which(xx=="G")
	res[k,se,4]=1	
}
dat.neg = res



#---------------------------
nn1 = nrow(dat.pos)
se1 = sample(1:nn1)[1:(nn1/2)]
nn2 = nrow(dat.neg)
se2 = sample(1:nn2)[1:(nn2/2)]

tmp1 = dat.pos[se1, ,]
tmp2 = dat.neg[se2, ,]
tr.x = abind(tmp1, tmp2, along=1)
dim(tr.x)
tr.y = c(rep(1, dim(tmp1)[1]), rep(0, dim(tmp2)[1]))
idx = sample(1:length(tr.y))
tr.x = tr.x[idx,,]
tr.y = tr.y[idx]

tmp1 = dat.pos[-se1, ,]
tmp2 = dat.neg[-se2, ,]
te.x = abind(tmp1, tmp2, along=1)
dim(tr.x)
te.y = c(rep(1, dim(tmp1)[1]), rep(0, dim(tmp2)[1]))
idx = sample(1:length(tr.y))
te.x = te.x[idx,,]
te.y = te.y[idx]

#---------------------------
#  tr.y <- to_categorical(tr.y, 2)				## skip this for two-class classification


#tr.x <- array_reshape(tr.x, c(dim(tr.x)[1], 1,100, 4))
#te.x <- array_reshape(te.x, c(dim(te.x)[1], 1, 100, 4))
dim(tr.x)
dim(te.x)

#adjust to 20000 1 100 4

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model <- keras_model_sequential() 
model %>% 
  layer_conv_1d(filter = 320, kernel_size = 10, activation = 'relu', input_shape = c(100, 4)) %>% 
  layer_dropout(rate = 0.25) %>%
  layer_max_pooling_1d(pool_size = 2) %>%
  layer_flatten() %>%
  layer_dense(2000,activation = 'relu') %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_dense(200,activation = 'relu') %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_dense(units = 20, activation = 'relu') %>%
  layer_dropout(rate = 0.25) %>%
#  layer_dense(units = 2, activation = 'softmax')		## change units=1 and activation="sigmoid" for two-class classification
  layer_dense(units = 1, activation = 'sigmoid')

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#max pooling
#kernel size
#epoach
summary(model)

#set up checkpoints
checkpoint_path <- "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/CNN/Model/random_model/cp.ckpt"

# Create checkpoint callback
cp_callback <- callback_model_checkpoint(
  filepath = checkpoint_path,
  save_weights_only = TRUE,
  save_best_only = TRUE,
  verbose = 1
)

model %>% compile(
#  loss = 'categorical_crossentropy',							## loss = "binary_crossentropy" for for two-class classification
  loss = 'binary_crossentropy',							
  optimizer = optimizer_rmsprop(lr = 0.0001, decay = 1e-6),
  metrics = c('accuracy')
)


history <- model %>% fit(
  tr.x, tr.y, 
  epochs = 50, batch_size = 128, 
  validation_split = 0.2,
  callbacks = list(cp_callback), # pass callback to training
)

myoutf ="/lorax/chenglab/yanding/ATAC_seq_integration/Figures/CNN/CNN_training_random_1d.pdf"
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
0.66

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/CNN/GM12878_random_CNN.Rda"
save.image(myoutf)
###########################################################################################################
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##########################################################################################################
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

#set up checkpoints
checkpoint_path <- "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/CNN/Model/random_model/cp2.ckpt"

# Create checkpoint callback
cp_callback <- callback_model_checkpoint(
  filepath = checkpoint_path,
  monitor = "val_accuracy",
  save_weights_only = TRUE,
  save_best_only = TRUE,
  verbose = 1
)

model %>% compile(
#  loss = 'categorical_crossentropy',							## loss = "binary_crossentropy" for for two-class classification
  loss = 'binary_crossentropy',							
  optimizer = optimizer_rmsprop(lr = 0.0001, decay = 1e-6),
  metrics = c('accuracy')
)


history <- model %>% fit(
  tr.x, tr.y, 
  epochs = 50, batch_size = 128, 
  validation_split = 0.2,
  callbacks = list(cp_callback), # pass callback to training
)


myoutf ="/lorax/chenglab/yanding/ATAC_seq_integration/Figures/CNN/CNN_training_random_1d.pdf"
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
0.66

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/CNN/GM12878_random_CNN_V2.Rda"
save.image(myoutf)
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[4]Model tuning 4th time
model <- keras_model_sequential() 
model %>% 
  layer_conv_1d(filter = 320, kernel_size = 5, activation = 'relu', input_shape = c(100, 4)) %>% 
  layer_conv_1d(filter = 320, kernel_size = 5, activation = 'relu', input_shape = c(100, 4)) %>% 
  layer_dropout(rate = 0.2) %>%
  layer_max_pooling_1d(pool_size = 2) %>%
  layer_conv_1d(filter = 480, kernel_size = 5, activation = 'relu', input_shape = c(100, 4)) %>% 
  layer_conv_1d(filter = 480, kernel_size = 5, activation = 'relu', input_shape = c(100, 4)) %>% 
  layer_dropout(rate = 0.2) %>%
  layer_max_pooling_1d(pool_size = 2) %>%
  layer_conv_1d(filter = 640, kernel_size = 5, activation = 'relu', input_shape = c(100, 4)) %>%  
  layer_flatten() %>%
  layer_dropout(rate = 0.2) %>% 
  layer_dense(2000,activation = 'relu') %>% 
  layer_dense(120,activation = 'relu') %>% 
#  layer_dense(units = 2, activation = 'softmax')		## change units=1 and activation="sigmoid" for two-class classification
  layer_dense(units = 1, activation = 'sigmoid')



#max pooling
#kernel size
#epoach
summary(model)

#set up checkpoints
checkpoint_path <- "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/CNN/Model/random_model/cp3.ckpt"

# Create checkpoint callback
cp_callback <- callback_model_checkpoint(
  filepath = checkpoint_path,
  monitor = "val_accuracy",
  save_weights_only = TRUE,
  save_best_only = TRUE,
  verbose = 1
)

model %>% compile(
#  loss = 'categorical_crossentropy',							## loss = "binary_crossentropy" for for two-class classification
  loss = 'binary_crossentropy',							
  optimizer = optimizer_rmsprop(lr = 0.0001, decay = 1e-6),
  metrics = c('accuracy')
)


history <- model %>% fit(
  tr.x, tr.y, 
  epochs = 50, batch_size = 128, 
  validation_split = 0.2,
  callbacks = list(cp_callback), # pass callback to training
)

myoutf ="/lorax/chenglab/yanding/ATAC_seq_integration/Figures/CNN/CNN_training_random_V3.pdf"
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
0.61

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/CNN/GM12878_random_CNN_V3.Rda"
save.image(myoutf)
