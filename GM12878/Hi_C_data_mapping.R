remotes::install_github("aidenlab/straw/R")

[1]Extract the contact information
rm(list=ls())
library(strawr)
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/Hi_C_hic/"
files = list.files(mydir)
sam = gsub(".hic","",files)

chr = c(seq(1,22, 1),"X")

for(i in 11 : length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir, files[i])
	
	for(k in 1 : length(chr))
	{
		cat("\r",i,"-->",k)
		hic.data.frame <- strawr::straw("VC", tmpinf, chr[k], chr[k], "BP", 50000)
		
		#delete the self loop
		tag = hic.data.frame$x == hic.data.frame$y
		hic.data.frame = hic.data.frame[tag==0,]
		
		hic.data.frame[,"chr"] = rep(paste0("chr",chr[k]),nrow(hic.data.frame))
		colnames(hic.data.frame) = c("sta_1","sta_2","counts","chr")
		hic.data.frame = hic.data.frame[,c("chr","sta_1","sta_2","counts")]
		hic.data.frame$end_1 = hic.data.frame$sta_1 + 10000
		hic.data.frame$end_2 = hic.data.frame$sta_2 + 10000
		hic.data.frame =  hic.data.frame[,c("chr","sta_1","end_1","sta_2","end_2","counts")]
		
		if(k == 1)
		{
			res = hic.data.frame
		}else{
			res = rbind(res, hic.data.frame)
		}
	}
	
	myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/Hi_C_contact_map/",sam[i],"_contact.txt")
	write.table(res,myoutf,sep="\t",quote=F)
}

[2]Begin to calculate the signal pos
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878_GM12878_Hi_C_pos/"
dir.create(myDir1)
mytemplate = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Template/Template_for_counting_GM12878_GM12878_Hi_C_sig_pos.R"
mysub = paste(myDir1, "submit.sp", sep="")
setwd(myDir1)

conIn = file(mytemplate, "r")
rawdata = readLines(conIn)
close(conIn)

fnum = 14
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

[5]Begin to calcualte the signal for negative
rm(list=ls())
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878_GM12878_Hi_C_neg/"
dir.create(myDir1)
mytemplate = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Template/Template_for_counting_GM12878_GM12878_Hi_C_sig_neg.R"
mysub = paste(myDir1, "submit.sp", sep="")
setwd(myDir1)

conIn = file(mytemplate, "r")
rawdata = readLines(conIn)
close(conIn)

fnum = 14
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