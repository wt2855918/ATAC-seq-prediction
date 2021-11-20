"/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_bin_sequence/ENCFF356TXH_positive_training.fa"
[1]Training data process
rm(list=ls())

#install.packages("tensorflow")
library(tensorflow)
#install_tensorflow()
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_bin_sequence/ENCFF356TXH_positive_training.fa"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_bin_sequence/ENCFF356TXH_negative_training.fa"

library(abind)
library(keras)
library("Biostrings")
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

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/CNN/HepG2_local_training_CNN.Rda"
save.image(file=myoutf)

##################################################################################
[2]Validation data process
rm(list=ls())

#install.packages("tensorflow")
library(tensorflow)
#install_tensorflow()
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_bin_sequence/ENCFF356TXH_positive_validation.fa"
myinf2 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_bin_sequence/ENCFF356TXH_negative_validation.fa"

library(abind)
library(keras)
library("Biostrings")
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
tmp1 = dat.pos
tmp2 = dat.neg
te.x = abind(tmp1, tmp2, along=1)
dim(tr.x)
te.y = c(rep(1, dim(tmp1)[1]), rep(0, dim(tmp2)[1]))

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/CNN/HepG2_local_validation_CNN.Rda"
save.image(file=myoutf)