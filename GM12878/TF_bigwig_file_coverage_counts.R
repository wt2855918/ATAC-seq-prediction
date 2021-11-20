[1]Begin to calculate the coverage 
rm(list=ls())
work_mydir  = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF_from_bigWig/"
dir.create(work_mydir)
setwd(work_mydir)

sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

myoutdir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/"
dir.create(myoutdir)

for(k in 1 : length(bam_files))
{
	cat("\r",k)
	myoutf1 = paste("job", k, ".sp", sep="")
	conOut = file(myoutf1, "w")
	curLine = c("#PBS -l walltime=10:00:00 -l vmem=256gb", 
				"", 
				"module load python/2.7-Anaconda", 
				"source activate SICER")
	writeLines(curLine, conOut)
	
	curLine = c("cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/")
	writeLines(curLine, conOut)
	
	input_file = bam_files[k]
	
	output1 = paste0(myoutdir,sam[k],".npz")
	output2 = paste0(myoutdir,sam[k],".tab")
	
	curLine = paste0("multiBigwigSummary bins -bs 100 -b ", input_file," -out ", output1, " --outRawCounts ",output2)
	
	writeLines(curLine, conOut)	
	close(conOut)
}

#begin to submit job
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF_from_bigWig/"
files = list.files(myDir1)

for(i in 1 : length(files))
{
	cat("\r",i)
	myinf = paste0(myDir1,files[i])
	command = paste0("qsub ", myinf)
	system(command)
}

##begin to check
#library(rtracklayer)
#myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/Histone_ChIP_seq_bigWig/ENCFF091LGA.bigWig"
#res = import(myinf, format = "bigWig", as = "GRanges")

[2]Examine the coverage that has not been submitted
##############################################
#Examine the bug files
##############################################
rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]
files = gsub(".tab","",files)

tag = which(sam %in% files)
com = sam[-tag]

tag = which(sam %in% com)
tar = paste0("job",tag,".sp")

work_mydir  = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF_from_bigWig/"
dir.create(work_mydir)
setwd(work_mydir)

#begin to submit job
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF_from_bigWig/"
files = list.files(myDir1)
files = tar
for(i in 1 : length(files))
{
	cat("\r",i)
	myinf = paste0(myDir1,files[i])
	command = paste0("qsub ", myinf)
	system(command)
}

[3]Some files are not still calculated even after increasing memory so decided to process with each individual chromsome
rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]
files = gsub(".tab","",files)

tag = which(sam %in% files)
com = sam[-tag]

for(i in 1 : length(com))
{
	cat("\r",i)
	input_file = paste0(com[i],".bigWig")
	
	scripts_dir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF_from_bigWig/",com[i])
	dir.create(scripts_dir)
	setwd(scripts_dir)

	chr_id = paste0("chr",c(seq(1,22,1),"X","Y"))
	
	myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i],"/")
	dir.create(myoutdir)
	
	
	for(k in 1 : length(chr_id))
	{
		myoutf1 = paste("job", k, ".sp", sep="")
		conOut = file(myoutf1, "w")
		curLine = c("#PBS -l walltime=10:00:00 -l vmem=256gb", 
				"", 
				"module load python/2.7-Anaconda", 
				"source activate SICER")
		writeLines(curLine, conOut)
		
		curLine = c("cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/")
		writeLines(curLine, conOut)
		
		output1 = paste0(myoutdir,com[i],"_",chr_id[k],".npz")
		output2 = paste0(myoutdir,com[i],"_",chr_id[k],".tab")
		curLine = paste0("multiBigwigSummary bins -bs 100 -r ",chr_id[k]," -b ", input_file," -out ", output1, " --outRawCounts ",output2)
		writeLines(curLine, conOut)	
		close(conOut)
	}
}

#begin to submit
for(i in 1: length(com))
{
	scripts_dir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF_from_bigWig/",com[i],"/")
	dir.create(scripts_dir)
	setwd(scripts_dir)
	
	files = list.files(scripts_dir)
	for(k in 1 : length(files))
	{
		cat("\r",i,"-->",k)
		myinf = paste0(scripts_dir,files[k])
		command = paste0("qsub ", myinf)
		system(command)
	}
}

[4]Some files that a single chromosome cannot be calculated
rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]
files = gsub(".tab","",files)

tag = which(sam %in% files)
com = sam[-tag]

chr.len = c(249250621, 243199373, 198022430, 191154276, 180915260,
            171115067, 159138663, 146364022, 141213431, 135534747,
            135006516, 133851895, 115169878, 107349540, 102531392,
            90354753, 90354753, 78077248, 59128983, 63025520,
            48129895, 51304566, 155270560,59373566)

chr_id = paste0("chr",c(seq(1,22,1),"X","Y"))
names(chr.len) = chr_id

for(i in 1 : length(com))
{
	cat("\r",i)
	tmpdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i])
	files = list.files(tmpdir)
	tag = grep(".tab",files)
	files = files[tag]
	nam = sapply(files,function(x) strsplit(x,"_")[[1]][2])
	nam = as.vector(gsub(".tab","",nam))
	
	tag = which(chr_id %in% nam)
	target = chr_id[-tag]
	
	input_file = paste0(com[i],".bigWig")
	
	if(length(target) > 0)
	{
		for(k in 1 : length(target))
		{
			scripts_dir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF_from_bigWig/",com[i],"/",target[k],"/")
			dir.create(scripts_dir)
			setwd(scripts_dir)
		
			myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i],"/",target[k],"/")
			dir.create(myoutdir)
		
			total_len = chr.len[target[k]]
			
			bin.size = 100000
			bin_num = ceiling(total_len/bin.size)
			
			for(j in 1 : bin_num)
			{
				cat("\r",i,"-->",k,"-->",j)
				sta = (j-1)*bin.size
				end = min(sta+bin.size,total_len)
				
				sta = as.integer(sta)
				end = as.integer(end)
				
				myoutf1 = paste("job", j, ".sp", sep="")
				conOut = file(myoutf1, "w")
				curLine = c("#PBS -l walltime=10:00:00 -l vmem=16gb", 
				"", 
				"module load python/2.7-Anaconda", 
				"source activate SICER")
				writeLines(curLine, conOut)
				
				curLine = c("cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/")
				writeLines(curLine, conOut)
				
				output1 = paste0(myoutdir,com[i],"_",target[k],"_",j,".npz")
				output2 = paste0(myoutdir,com[i],"_",target[k],"_",j,".tab")
				
				curLine = paste0("multiBigwigSummary bins -bs 100 -r ",target[k],":",sta,":",end," -b ", input_file," -out ", output1, " --outRawCounts ",output2)
				writeLines(curLine, conOut)	
				close(conOut)
			
			}
		
		}
	}else{
		next
	}

}

[5]Even at the 10k level, some files are failed...examine again..
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF_from_bigWig/ENCFF285JUZ/chr8/"

rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]
files = gsub(".tab","",files)

tag = which(sam %in% files)
com = sam[-tag]

chr.len = c(249250621, 243199373, 198022430, 191154276, 180915260,
            171115067, 159138663, 146364022, 141213431, 135534747,
            135006516, 133851895, 115169878, 107349540, 102531392,
            90354753, 90354753, 78077248, 59128983, 63025520,
            48129895, 51304566, 155270560,59373566)

chr_id = paste0("chr",c(seq(1,22,1),"X","Y"))
names(chr.len) = chr_id

for(i in 1 : length(com))
{
	cat("\r",i)
	tmpdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i])
	files = list.files(tmpdir)
	tag = grep(".tab",files)
	files = files[tag]
	nam = sapply(files,function(x) strsplit(x,"_")[[1]][2])
	nam = as.vector(gsub(".tab","",nam))
	
	tag = which(chr_id %in% nam)
	target = chr_id[-tag]
	
	input_file = paste0(com[i],".bigWig")
	
	if(length(target) > 0)
	{
		for(k in 1 : length(target))
		{
			scripts_dir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF_from_bigWig/",com[i],"/",target[k],"/")
			dir.create(scripts_dir)
			setwd(scripts_dir)
		
			myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i],"/",target[k],"/")
			dir.create(myoutdir)
		
			total_len = chr.len[target[k]]
			
			bin.size = 100000
			bin_num = ceiling(total_len/bin.size)
			
			sub_job_ID = paste0("job",seq(1,bin_num,1),".sp")
			
			job_finish = list.files(myoutdir)
			tag = grep(".tab",job_finish)
			job_finish = job_finish[tag]
			job_finish = sapply(job_finish, function(x) strsplit(x,"_")[[1]][3])
			job_finish = gsub(".tab","",job_finish)
			job_finish = paste0("job",job_finish,".sp")
			
			tag = which(sub_job_ID %in% job_finish)
			sub_job_ID = sub_job_ID[-tag]
			
			for(j in 1 : length(sub_job_ID))
			{
				cat("\r",i,"-->",k,"-->",j)
				
				conIn = file(sub_job_ID[j], "r")
				rawdata = readLines(conIn)
				close(conIn)
				
				rawdata[1] = c("#PBS -l walltime=10:00:00 -l vmem=128gb")
				
				myoutf1 = sub_job_ID[j]
				conOut = file(myoutf1, "w")
				writeLines(rawdata, conOut)
				close(conOut)
				
				command = paste0("qsub ",sub_job_ID[j])
				system(command)
			}
		
		}
	}else{
		next
	}

}

[6]still have issues ...go to 10000 base pair layer
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF_from_bigWig/ENCFF285JUZ/chr8/"

rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]
files = gsub(".tab","",files)

tag = which(sam %in% files)
com = sam[-tag]

chr.len = c(249250621, 243199373, 198022430, 191154276, 180915260,
            171115067, 159138663, 146364022, 141213431, 135534747,
            135006516, 133851895, 115169878, 107349540, 102531392,
            90354753, 90354753, 78077248, 59128983, 63025520,
            48129895, 51304566, 155270560,59373566)

chr_id = paste0("chr",c(seq(1,22,1),"X","Y"))
names(chr.len) = chr_id

for(i in 1 : length(com))
{
	tmpdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i])
	files = list.files(tmpdir)
	tag = grep(".tab",files)
	files = files[tag]
	nam = sapply(files,function(x) strsplit(x,"_")[[1]][2])
	nam = as.vector(gsub(".tab","",nam))
	
	tag = which(chr_id %in% nam)
	target = chr_id[-tag]
	
	input_file = paste0(com[i],".bigWig")
	
	if(length(target) > 0)
	{
		for(k in 1 : length(target))
		{
			scripts_dir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF_from_bigWig/",com[i],"/",target[k],"/")
			dir.create(scripts_dir)
			setwd(scripts_dir)
		
			myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i],"/",target[k],"/")
			dir.create(myoutdir)
		
			total_len = chr.len[target[k]]
			
			bin.size = 100000
			bin_num = ceiling(total_len/bin.size)
			
			sub_job_ID = paste0("job",seq(1,bin_num,1),".sp")
			
			job_finish = list.files(myoutdir)
			tag = grep(".tab",job_finish)
			job_finish = job_finish[tag]
			job_finish = sapply(job_finish, function(x) strsplit(x,"_")[[1]][3])
			job_finish = gsub(".tab","",job_finish)
			job_finish = paste0("job",job_finish,".sp")
			
			tag = which(sub_job_ID %in% job_finish)
			sub_job_ID = sub_job_ID[-tag]
			
			for(j in 1 : length(sub_job_ID))
			{
			
				scripts_dir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF_from_bigWig/",com[i],"/",target[k],"/")
				dir.create(scripts_dir)
				setwd(scripts_dir)
		
				myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i],"/",target[k],"/")
				dir.create(myoutdir)
		
				conIn = file(sub_job_ID[j], "r")
				rawdata = readLines(conIn)
				close(conIn)
				
				position_line = rawdata[6]
				position_line = strsplit(position_line," ")[[1]][6]
				
				sta_p = strsplit(position_line,":")[[1]][2]
				end_p = strsplit(position_line,":")[[1]][3]
				
				sta_p = as.numeric(sta_p)
				end_p = as.numeric(end_p)
				
				total_len = as.numeric(end_p) - as.numeric(sta_p)
				bin.size = 10000
				
				bin_num = ceiling(total_len/bin.size)
				
				for(m in 1 : bin_num)
				{
					cat("\r",i,"-->",k,"-->",j,"-->",m)
					sta = (m-1)*bin.size + sta_p
					end = min(sta+bin.size,end_p)
					
					sta = as.integer(sta)
					end = as.integer(end)
					
					myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i],"/",target[k],"/")
	
					myoutdir = paste0(myoutdir,"sub_bin/")
					dir.create(myoutdir)
						
					output1 = paste0(myoutdir,com[i],"_",target[k],"_",j,".npz")
					output2 = paste0(myoutdir,com[i],"_",target[k],"_",j,".tab")
				
					curLine =  paste0("multiBigwigSummary bins -bs 100 -r ",target[k],":",sta,":",end," -b ", input_file," -out ", output1, " --outRawCounts ",output2)
					
					rawdata[6] = curLine
					rawdata[1] = c("#PBS -l walltime=10:00:00 -l vmem=256gb")
					
					scripts_dir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878/Counting_TF_from_bigWig/",com[i],"/",target[k],"/sub_bin/")
					dir.create(scripts_dir)
					setwd(scripts_dir)
					
					old_job_tag = gsub(".sp","",sub_job_ID[j])
					
					myoutf1 = paste0(old_job_tag,"_",m,".sp")
					conOut = file(myoutf1, "w")
					writeLines(rawdata, conOut)
					close(conOut)
				}

			}
		
		}
	}else{
		next
	}

}

module load python/2.7-Anaconda
source activate SICER
cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/
multiBigwigSummary bins -bs 100 -r chr8 -b ENCFF285JUZ.bigWig -out /lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/ENCFF285JUZ/ENCFF285JUZ_chr8.npz --outRawCounts /lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/ENCFF285JUZ/ENCFF285JUZ_chr8.tab

#mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/ENCFF942IMP/"
#myoutdir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/ENCFF942IMP/chr8/"

#files = list.files(mydir)
#tag = grep("chr8",files)
#files = files[tag]
#tag = grep(".tab",files)
#files = files[tag]

#for(i in 1 : length(files))
#{
#	cat("\r",i)
#	inputfiles = paste0(mydir,files[i])
#	outputfiles = paste0(myoutdir,files[i])
#	command = paste0("mv ",inputfiles," ",outputfiles)
#	system(command)
#}

[6]Begin to merge a single chromosome subbin (10k level)
rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]
files = gsub(".tab","",files)

tag = which(sam %in% files)
com = sam[-tag]

chr.len = c(249250621, 243199373, 198022430, 191154276, 180915260,
            171115067, 159138663, 146364022, 141213431, 135534747,
            135006516, 133851895, 115169878, 107349540, 102531392,
            90354753, 90354753, 78077248, 59128983, 63025520,
            48129895, 51304566, 155270560,59373566)

chr_id = paste0("chr",c(seq(1,22,1),"X","Y"))
names(chr.len) = chr_id

for(i in 1 : length(com))
{
	tmpdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i])
	files = list.files(tmpdir)
	tag = grep(".tab",files)
	files = files[tag]
	nam = sapply(files,function(x) strsplit(x,"_")[[1]][2])
	nam = as.vector(gsub(".tab","",nam))
	
	tag = which(chr_id %in% nam)
	target = chr_id[-tag]
	
	input_file = paste0(com[i],".bigWig")
	
	if(length(target) > 0)
	{
		for(k in 1 : length(target))
		{
			tmpinf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i],"/",target[k],"/sub_bin/")
			bin_files = list.files(tmpinf)
			tag = grep(".tab",bin_files)
			bin_files = bin_files[tag]
			
			myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i],"/",target[k],"/")
			
			for(j in 1 : length(bin_files))
			{
				input_files = paste0(tmpinf,bin_files[j])
				data = read.table(input_files,sep="\t",quote=NULL)
				
				
				my_nam = gsub(".tab","",bin_files[j])
				my_nam = paste0(my_nam,"_",j,".tab")
				
				output_files = paste0(myoutdir,my_nam)
				
				command = paste0("mv ",input_files," ",output_files)
				system(command)
			}
		}
	}
}	

[7]Summary a chromosome
rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]
files = gsub(".tab","",files)

tag = which(sam %in% files)
com = sam[-tag]

chr.len = c(249250621, 243199373, 198022430, 191154276, 180915260,
            171115067, 159138663, 146364022, 141213431, 135534747,
            135006516, 133851895, 115169878, 107349540, 102531392,
            90354753, 90354753, 78077248, 59128983, 63025520,
            48129895, 51304566, 155270560,59373566)

chr_id = paste0("chr",c(seq(1,22,1),"X","Y"))
names(chr.len) = chr_id

for(i in 1 : length(com))
{
	tmpdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i])
	files = list.files(tmpdir)
	tag = grep(".tab",files)
	files = files[tag]
	nam = sapply(files,function(x) strsplit(x,"_")[[1]][2])
	nam = as.vector(gsub(".tab","",nam))
	
	tag = which(chr_id %in% nam)
	target = chr_id[-tag]
	
	input_file = paste0(com[i],".bigWig")
	
	if(length(target) > 0)
	{
		for(k in 1 : length(target))
		{
			tmpinf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i],"/",target[k],"/")
			bin_files = list.files(tmpinf)
			tag = grep(".tab",bin_files)
			bin_files = bin_files[tag]
			
			myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i],"/")
			
			for(j in 1 : length(bin_files))
			{
				cat("\r",i,"-->",k,"-->",j)
				input_files = paste0(tmpinf,bin_files[j])
				data = read.table(input_files,sep="\t",quote=NULL)
				
				if(j == 1)
				{
				
					res = data
				}else{
				
					res = rbind(data,res)
				
				}
				
				myoutf = paste0(myoutdir,"/",com[i],"_",target[k],".tab")
				write.table(res,myoutf,sep="\t",quote=F)
				
			}
		}
	}
}	


[8]Begin to merge into a single file
rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]
files = gsub(".tab","",files)

tag = which(sam %in% files)
com = sam[-tag]

for(i in 1 : length(com))
{
	myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i],".tab")

	chr_id = paste0("chr",c(seq(1,22,1),"X","Y"))
	
	for(k in 1 : length(chr_id))
	{
		cat("\r",i,"-->",k)
		tmpinf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/",com[i],"/",com[i],"_",chr_id[k],".tab")
		res = read.table(tmpinf,sep="\t",quote=NULL)
		add = c(chr_id[k],"0","100","0")
		res = rbind(add,res)
		
		if(k==1)
		{
		
			final_res = res
		}else{
		
			final_res = rbind(final_res,res)
		}
	}

	write.table(final_res,myoutf,sep="\t",quote=F)
}

[5]Begin to calculate the signal
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]

mytardir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_counts_bins_logFC/"
my_end = list.files(mytardir)
tag = grep("positive",my_end)
my_end = my_end[tag]
my_end = gsub("_ENCFF172DEA_positive.Rda","",my_end)
my_end = unique(my_end)

my_nam = gsub(".tab","",files)
tag = which(my_nam %in% my_end)
if(length(tag)>0)
{
	files = files[-tag]
}else{
	files = files
}
fnum = length(files)

myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878_GM12878_TF_pos/"
dir.create(myDir1)
mytemplate = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Template/Template_for_counting_GM12878_GM12878_TF_sig_pos.R"
mysub = paste(myDir1, "submit.sp", sep="")
setwd(myDir1)

conIn = file(mytemplate, "r")
rawdata = readLines(conIn)
close(conIn)

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
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]

mytardir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_counts_bins_logFC/"
my_end = list.files(mytardir)
tag = grep("negative",my_end)
my_end = my_end[tag]
my_end = gsub("_ENCFF172DEA_negative.Rda","",my_end)
my_end = unique(my_end)

my_nam = gsub(".tab","",files)
tag = which(my_nam %in% my_end)
if(length(tag)>0)
{
	files = files[-tag]
}else{
	files = files
}
fnum = length(files)

myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/GM12878_GM12878_TF_neg/"
dir.create(myDir1)
mytemplate = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Template/Template_for_counting_GM12878_GM12878_TF_sig_neg.R"
mysub = paste(myDir1, "submit.sp", sep="")
setwd(myDir1)

conIn = file(mytemplate, "r")
rawdata = readLines(conIn)
close(conIn)

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

[5]Model fiting
#begin to test the shared number (Positive)
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_counts_bins_logFC/"
files = list.files(mydir)
tag = grep("positive",files)
files = files[tag]

for(i in 1: length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir, files[i])
	load(tmpinf)
	
	if(i == 1)
	{
		bins = histone$ID
	}else{
	
		bins = intersect(bins, histone$ID)
	}
	
	if(length(bins) < 1)
	{
		command = paste0("Warning job",i)
		print(command)
	}
}
length(bins)
[1] 895695

#begin to test the shared number (Positive)
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_counts_bins_logFC/"
files = list.files(mydir)
tag = grep("negative",files)
files = files[tag]

for(i in 1: length(files))
{
	cat("\r",i)
	tmpinf = paste0(mydir, files[i])
	load(tmpinf)
	
	if(i == 1)
	{
		bins = histone$ID
	}else{
	
		bins = intersect(bins, histone$ID)
	}

}
length(bins)
[1] 895704

rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_counts_bins_logFC/"
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/TF_file_annotation_bigWig.txt"

files = list.files(mydir)
name = gsub(".Rda","",files)
name = gsub("_ENCFF172DEA_positive","",name)
name = gsub("_ENCFF172DEA_negative","",name)
name = unique(name)

info = read.table(myinf1,sep="\t",quote=NULL,stringsAsFactors=F)

tag = which(row.names(info) %in% name)
info = info[intersect(name,row.names(info)),]

library(randomForest)
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Multi_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Uni_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/moveme.R")

AUC_score = list()

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_positive.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_negative.Rda")


tag_pos = sample(seq(1,895695),2000)
tag_neg = sample(seq(1,895704),2000)

rm(positive)
rm(negative)

res = matrix(0, 4000, nrow(info))
colnames(res) = row.names(info)
res = as.data.frame(res)

for(i in 1 : nrow(info))
{
	cat("\r",i)
	
	tmpinf1 = paste0(mydir,row.names(info)[i],"_ENCFF172DEA_negative.Rda")
	load(tmpinf1)
	negative = histone
	negative = negative[tag_neg,]
	
	tmpinf2 = paste0(mydir,row.names(info)[i],"_ENCFF172DEA_positive.Rda")
	load(tmpinf2)
	positive = histone
	positive = positive[tag_pos,]
	
	negative$tag = rep(0,nrow(negative))
	positive$tag = rep(1,nrow(positive))
	
	data = rbind(positive,negative)
	
	if(i == 1)
	{
		row.names(res) = row.names(data)
		res[,i] = data$V4
	}else{
	
		res[,i] = data$V4
	}

}
raw.res = res
res = apply(res,2,function(x) as.numeric(as.vector(x)))
res = sign(res)
tag = c(rep(1,2000),rep(0,2000))
data = cbind(tag,res)
data = as.data.frame(data)
data$tag = factor(data$tag, levels =c(0,1))

source("/lorax/chenglab/yanding/ATAC_seq_integration/Function/RF_multi_var_prediction_10_CV.R")
fit = ROC(data)

1-fit[[3]]
[1] 0.8772474 0.8775325 0.8769620 0.8775456 0.8768321 0.8778768 0.8777907
 [8] 0.8767283 0.8765690 0.8774515
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_TF_integration_AUC.Rda"
save.image(file = myoutf)

