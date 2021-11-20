[2]Begin to calculate shared bins template code
rm(list=ls())
#define the function
f1 <- function(x, chr, sta, end) {paste(x[chr], x[sta]:x[end], sep="_")} #separate each base pair
f2 <- function(x, sig, sta, end) 
	{
		time = x[end]-x[sta]+1
		rep(as.character(x[sig]), time)
	}
	
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Histone_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]


myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/ATAC_seq_peak_processed_bed/ENCFF356TXH_positive.Rda"
load(myinf)
bin.size = 100

mytardir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Histone_counts_bins_logFC_GM12878_cross_HepG2/"
my_nam = gsub(".tab","",files)

i= mykkk
cat("\r",i)
tmpinf = paste0(mydir,files[i])
data = read.table(tmpinf,sep="\t",quote=NULL)
data[,"V3"] = data[,"V3"] - 1

#delete the uncoverage region (due to hg19 version difference)
tag = is.na(data[,"V4"])
data = data[tag==0,]
	
tar = paste0("chr",c(seq(1,22,1),"X","Y"))
tag = which(data[,"V1"] %in% tar)
data = data[tag,]
	
histone = data
histone[,"ID"] = paste0(histone[,"V1"],"_",histone[,"V2"],"_",histone[,"V3"])

positive[,"ID"] = paste0(positive[,"V1"],"_",positive[,"V2"],"_",positive[,"V3"])

row.names(histone) = histone[,"ID"]
row.names(positive) = positive[,"ID"]
	
com = intersect(row.names(histone),row.names(positive))
	
histone = histone[com,]
	
myoutf = paste0(mytardir,my_nam[i],"_ENCFF172DEA_positive.Rda")
save(histone,file = myoutf)

[3]Begin to calculate the signal
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878_HepG2_Histone/"
dir.create(myDir1)
mytemplate = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Template/Template_for_counting_GM12878_GM12878_Histone_sig.R"
mysub = paste(myDir1, "submit.sp", sep="")
setwd(myDir1)

conIn = file(mytemplate, "r")
rawdata = readLines(conIn)
close(conIn)

fnum = 39
for(k in 1:fnum)
{
	data = rawdata
	se = grep("i= mykkk", data)
	tmp = data[se]
	tmp = gsub("mykkk", k, tmp)
	data[se] = tmp
	myoutf1 = paste("job", k, ".sp", sep="")
	conOut = file(myoutf1, "w")
	writeLines(data, conOut)
	close(conOut)
}

[4]Begin to calcualte the signal for negative
rm(list=ls())
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878_GM12878_Histone/"
dir.create(myDir1)
mytemplate = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Template/Template_for_couting_GM12878_GM12878_Histone_sig_neg.R"
mysub = paste(myDir1, "submit.sp", sep="")
setwd(myDir1)

conIn = file(mytemplate, "r")
rawdata = readLines(conIn)
close(conIn)

fnum = 39
for(k in 1:fnum)
{
	data = rawdata
	se = grep("i= mykkk", data)
	tmp = data[se]
	tmp = gsub("mykkk", k, tmp)
	data[se] = tmp
	myoutf1 = paste("job", k, ".sp", sep="")
	conOut = file(myoutf1, "w")
	writeLines(data, conOut)
	close(conOut)
}
