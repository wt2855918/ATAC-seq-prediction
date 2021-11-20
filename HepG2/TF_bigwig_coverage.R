[1]Begin to calculate the coverage 
rm(list=ls())
work_mydir  = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/"
dir.create(work_mydir)
setwd(work_mydir)

sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

myoutdir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
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
	
	curLine = c("cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/")
	writeLines(curLine, conOut)
	
	input_file = bam_files[k]
	
	output1 = paste0(myoutdir,sam[k],".npz")
	output2 = paste0(myoutdir,sam[k],".tab")
	
	curLine = paste0("multiBigwigSummary bins -bs 100 -b ", input_file," -out ", output1, " --outRawCounts ",output2)
	
	writeLines(curLine, conOut)	
	close(conOut)
}

#begin to submit job
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/"
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
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]
files = gsub(".tab","",files)

tag = which(sam %in% files)
com = sam[-tag]

tag = which(sam %in% com)
tar = paste0("job",tag,".sp")

work_mydir  = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/"
dir.create(work_mydir)
setwd(work_mydir)

#begin to submit job
myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/"
files = list.files(myDir1)
files = tar
for(i in 1 : length(tar))
{
	cat("\r",i)
	myinf = paste0(myDir1,tar[i])
	command = paste0("qsub ", myinf)
	system(command)
}
tar
[1] "job141.sp" "job437.sp"

[1]"ENCFF580GNN" "ENCFF938OUY" "ENCFF993RSW"


#covert to bedgraph
cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/
/lorax/chenglab/yanding/Software/bigWigToBedGraph ENCFF580GNN.bigWig /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/test/ENCFF580GNN.bedgraph
/lorax/chenglab/yanding/Software/bigWigToBedGraph ENCFF938OUY.bigWig /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/test/ENCFF938OUY.bedgraph
/lorax/chenglab/yanding/Software/bigWigToBedGraph ENCFF993RSW.bigWig /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/test/ENCFF993RSW.bedgraph

#clean bedgraph
rm(list=ls())
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/test/"
files = list.files(myinf)
tag = grep(".bedgraph",files)
files = files[tag]
tag = grep("_clean",files)
files = files[-tag]

my_nam = gsub(".bedgraph","",files)

for(i in 1 : length(files))
{
	cat("\r",i)
	
	tmpinf = paste0(myinf,files[i])
	info = read.table(tmpinf,sep="\t",quote=NULL)
	
	tar = paste0("chr",c(seq(1,22,1),"X","Y"))
	tag = which(info$V1 %in% tar)
	info = info[tag,]

	myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/test/",my_nam[i],"_clean.bedgraph")
	write.table(info,myoutf,sep="\t",quote=F,row.names=F,col.names=F)
}

cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/test/
/lorax/chenglab/yanding/Software/bedGraphToBigWig ENCFF580GNN_clean.bedgraph hg19_chrom_length.txt ENCFF580GNN.bigWig
/lorax/chenglab/yanding/Software/bedGraphToBigWig ENCFF938OUY_clean.bedgraph hg19_chrom_length.txt ENCFF938OUY.bigWig
/lorax/chenglab/yanding/Software/bedGraphToBigWig ENCFF993RSW_clean.bedgraph hg19_chrom_length.txt ENCFF993RSW.bigWig


[3]Some files are not still calculated even after increasing memory so decided to process with each individual chromsome
rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)


mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
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
	
	scripts_dir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/",com[i])
	dir.create(scripts_dir)
	setwd(scripts_dir)

	chr_id = paste0("chr",c(seq(1,22,1),"X","Y"))
	
	myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i],"/")
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
		
		curLine = c("cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/")
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
	scripts_dir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/",com[i],"/")
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

#All files have chromosome 8 issues,lets calcualte using the new files
cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/test/

module load python/2.7-Anaconda
source activate SICER
multiBigwigSummary bins -bs 100 -r chr8 -b ENCFF580GNN.bigWig -out /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/ENCFF580GNN/ENCFF580GNN_chr8.npz --outRawCounts /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/ENCFF580GNN/ENCFF580GNN_chr8.tab
multiBigwigSummary bins -bs 100 -r chr8 -b ENCFF938OUY.bigWig -out /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/ENCFF938OUY/ENCFF938OUY_chr8.npz --outRawCounts /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/ENCFF938OUY/ENCFF938OUY_chr8.tab
multiBigwigSummary bins -bs 100 -r chr8 -b ENCFF993RSW.bigWig -out /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/ENCFF993RSW/ENCFF993RSW_chr8.npz --outRawCounts /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/ENCFF993RSW/ENCFF993RSW_chr8.tab


[4]Some chromsomes are not calculated use 100000 base pair instead
rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
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
	tmpdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i])
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
			scripts_dir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/",com[i],"/",target[k],"/")
			dir.create(scripts_dir)
			setwd(scripts_dir)
		
			myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i],"/",target[k],"/")
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
				
				curLine = c("cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/")
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

#examine the different folders
#[1]/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF015ETD/chr8/ #done
#[2]/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF387OZD/chr8/ #done
#[3]/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF541SRY/chr8/ #done
#[4]/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF580GNN/chr8/ #done
#[5]/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF833UKH/chr8/ #done
#[6]/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF938OUY/chr8/ #done
#[7]/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF993RSW/chr8/ #done
####


[5]still have issues ...go to 10000 base pair layer

rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
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
	tmpdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i])
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
			scripts_dir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/",com[i],"/",target[k],"/")
			dir.create(scripts_dir)
			setwd(scripts_dir)
		
			myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i],"/",target[k],"/")
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
			
				scripts_dir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/",com[i],"/",target[k],"/")
				dir.create(scripts_dir)
				setwd(scripts_dir)
		
				myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i],"/",target[k],"/")
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
					
					myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i],"/",target[k],"/")
	
					myoutdir = paste0(myoutdir,"sub_bin/")
					dir.create(myoutdir)
						
					output1 = paste0(myoutdir,com[i],"_",target[k],"_",j,".npz")
					output2 = paste0(myoutdir,com[i],"_",target[k],"_",j,".tab")
				
					curLine =  paste0("multiBigwigSummary bins -bs 100 -r ",target[k],":",sta,":",end," -b ", input_file," -out ", output1, " --outRawCounts ",output2)
					
					rawdata[6] = curLine
					rawdata[1] = c("#PBS -l walltime=10:00:00 -l vmem=16gb")
					
					scripts_dir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/",com[i],"/",target[k],"/sub_bin/")
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



#begin to submit
rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
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
	tmpdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i])
	files = list.files(tmpdir)
	tag = grep(".tab",files)
	files = files[tag]
	nam = sapply(files,function(x) strsplit(x,"_")[[1]][2])
	nam = as.vector(gsub(".tab","",nam))
	
	tag = which(chr_id %in% nam)
	target = chr_id[-tag]
	
	input_file = paste0(com[i],".bigWig")
	
	for(k in 1 : length(target))
	{
		scripts_dir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/",com[i],"/",target[k],"/sub_bin/")
		files = list.files(scripts_dir)
		setwd(scripts_dir)
		
		for(j in 1 : length(files))
		{
			cat("\r",i,"-->",k,"-->",j)
			command = paste0("qsub ",files[j])
			system(command)
		
		}
	}
	
}


#PBS -l walltime=10:00:00 -l vmem=16gb

module load python/2.7-Anaconda
source activate SICER
cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/
multiBigwigSummary bins -bs 100 -r chr8:14580000:145810000 -b ENCFF015ETD.bigWig -out /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/ENCFF015ETD/chr8/sub_bin/ENCFF015ETD_chr8_10.npz --outRawCounts /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/ENCFF015ETD/chr8/sub_bin/ENCFF015ETD_chr8_10.tab
job1460_3.sp (END)

min(146364022, 145800000)


#examine the different folders
cd /lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF015ETD/chr8/sub_bin/ #done
rm -f *
cd /lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF387OZD/chr8/sub_bin/ #done
rm -f *
cd /lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF541SRY/chr8/sub_bin/ #done
rm -f *
cd /lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF580GNN/chr8/sub_bin/ #done
rm -f *
cd /lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF833UKH/chr8/sub_bin/ #done
rm -f *
cd /lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF938OUY/chr8/sub_bin/ #done
rm -f *
cd /lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF993RSW/chr8/sub_bin/ #done
rm -f *

cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/ENCFF015ETD/chr8/sub_bin/

####
#examine the different folders
#[1]/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF015ETD/chr8/ #done
#[2]/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF387OZD/chr8/ #done
#[3]/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF541SRY/chr8/ #done
#[4]/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF580GNN/chr8/ #done
#[5]/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF833UKH/chr8/ #done
#[6]/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF938OUY/chr8/ #done
#[7]/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/Counting_TF_from_bigWig/ENCFF993RSW/chr8/ #done
####
cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/
rm -f ENCFF015ETD.tab
rm -f ENCFF387OZD.tab
rm -f ENCFF541SRY.tab
rm -f ENCFF580GNN.tab
rm -f ENCFF833UKH.tab
rm -f ENCFF938OUY.tab
rm -f ENCFF993RSW.tab

[6]Begin to merge a single chromosome subbin (10k level)
rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
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
	tmpdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i])
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
			tmpinf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i],"/",target[k],"/sub_bin/")
			bin_files = list.files(tmpinf)
			tag = grep(".tab",bin_files)
			bin_files = bin_files[tag]
			myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i],"/",target[k],"/")
			
			if(length(bin_files)>0)
			{
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
}	

[7]Summary a chromosome
rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
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
	tmpdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i])
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
			tmpinf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i],"/",target[k],"/")
			bin_files = list.files(tmpinf)
			tag = grep(".tab",bin_files)
			bin_files = bin_files[tag]
			
			myoutdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i],"/")
			
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

[8]Summary a single file
rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]
files = gsub(".tab","",files)

tag = which(sam %in% files)
com = sam[-tag]

chr_id = paste0("chr",c(seq(1,22,1),"X","Y"))
names(chr.len) = chr_id

for(i in 1 : length(com))
{
	tmpdir = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i],"/")
	files = list.files(tmpdir)
	tag = grep(".tab",files)
	files = files[tag]

	for(k in 1 : length(files))
	{
		cat("\r",i,"-->",k)
		tmpinf = paste0(tmpdir, files[k])
		data = read.table(tmpinf,sep="\t",quote=NULL)
		
		if(k == 1)
		{
			res = data
		
		}else{
		
			res = rbind(data,res)
		
		}
		
		myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i],".tab")
		write.table(res,myoutf,sep="\t",quote=F)
	}
}

[9]Begin to merge into a single file
rm(list=ls())
sourcedir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/"
bam_files = list.files(sourcedir)
tag = grep(".bigWig",bam_files)
bam_files = bam_files[tag]
sam = gsub(".bigWig","",bam_files)

mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]
files = gsub(".tab","",files)

tag = which(sam %in% files)
com = sam[-tag]

for(i in 1 : length(com))
{
	myoutf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i],".tab")

	chr_id = paste0("chr",c(seq(1,22,1),"X","Y"))
	
	for(k in 1 : length(chr_id))
	{
		cat("\r",i,"-->",k)
		tmpinf = paste0("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/",com[i],"/",com[i],"_",chr_id[k],".tab")
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

#process bugs
cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/
rm -f ENCFF015ETD.tab
rm -f ENCFF387OZD.tab
rm -f ENCFF541SRY.tab
rm -f ENCFF580GNN.tab
rm -f ENCFF833UKH.tab
rm -f ENCFF938OUY.tab
rm -f ENCFF993RSW.tab

############################
myinf1 = "ENCFF015ETD.tab"
data = read.table(myinf1,sep="\t",quote=NULL)
tag = data$V1 == "chr8" & data$V2 == "0"
test = data[tag,]
data = data[-13913507,]
write.table(data,myinf1,sep="\t",quote=F)
############################
rm(list=ls())
myinf1 = "ENCFF387OZD.tab"
data = read.table(myinf1,sep="\t",quote=NULL)
tag = data$V1 == "chr8" & data$V2 == "0"
test = data[tag,]
xx = row.names(test)[1]
data = data[-13913507,]
write.table(data,myinf1,sep="\t",quote=F)
############################
rm(list=ls())
myinf1 = "ENCFF541SRY.tab"
data = read.table(myinf1,sep="\t",quote=NULL)
tag = data$V1 == "chr8" & data$V2 == "0"
test = data[tag,]
data = data[-13913507,]
write.table(data,myinf1,sep="\t",quote=F)
############################
rm(list=ls())
myinf1 = "ENCFF580GNN.tab"
data = read.table(myinf1,sep="\t",quote=NULL)
tag = data$V1 == "chr8" & data$V2 == "0"
test = data[tag,]
data = data[-13927961,]
write.table(data,myinf1,sep="\t",quote=F)
############################
rm(list=ls())
myinf1 = "ENCFF833UKH.tab"
data = read.table(myinf1,sep="\t",quote=NULL)
tag = data$V1 == "chr8" & data$V2 == "0"
test = data[tag,]
xx = row.names(test)[1]
data = data[-13913507,]
write.table(data,myinf1,sep="\t",quote=F)
############################
rm(list=ls())
myinf1 = "ENCFF938OUY.tab"
data = read.table(myinf1,sep="\t",quote=NULL)
tag = data$V1 == "chr8" & data$V2 == "0"
test = data[tag,]
xx = row.names(test)[1]
data = data[-13927961,]
write.table(data,myinf1,sep="\t",quote=F)
############################
rm(list=ls())
myinf1 = "ENCFF993RSW.tab"
data = read.table(myinf1,sep="\t",quote=NULL)
tag = data$V1 == "chr8" & data$V2 == "0"
test = data[tag,]
xx = row.names(test)[1]
data = data[-13927961,]
write.table(data,myinf1,sep="\t",quote=F)


[5]Begin to calculate the signal for pos
rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]

mytardir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_counts_bins_logFC/"
my_end = list.files(mytardir)
tag = grep("positive",my_end)
my_end = my_end[tag]
my_end = gsub("_ENCFF356TXH_positive.Rda","",my_end)
my_end = unique(my_end)

my_nam = gsub(".tab","",files)
tag = which(my_nam %in% my_end)

#if(length(tag)>0)
#{
#	files = files[-tag]
#}else{
#	files = files
#}

#induce the annotation file 
annotation_file_path = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/TF_file_annotation_bigWig.txt"
annotation = read.table(annotation_file_path,sep="\t",quote=NULL)
annotation_files = paste0(row.names(annotation),".tab")
files = intersect(annotation_files,files)

fnum = length(files)

myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2/HepG2_TF_pos/"
dir.create(myDir1)
mytemplate = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Template/HepG2/Template_for_counting_HepG2_TF_sig_pos.R"
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
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_logFC/"
files = list.files(mydir)
tag = grep(".tab",files)
files = files[tag]

mytardir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/TF_counts_bins_logFC/"
my_end = list.files(mytardir)
tag = grep("negative",my_end)
my_end = my_end[tag]
my_end = gsub("_ENCFF356TXH_negative.Rda","",my_end)
my_end = unique(my_end)

my_nam = gsub(".tab","",files)
tag = which(my_nam %in% my_end)

#if(length(tag)>0)
#{
#	files = files[-tag]
#}else{
#	files = files
#}

#induce the annotation file 
annotation_file_path = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/TF_file_annotation_bigWig.txt"
annotation = read.table(annotation_file_path,sep="\t",quote=NULL)
annotation_files = paste0(row.names(annotation),".tab")
files = intersect(annotation_files,files)

fnum = length(files)

myDir1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/HepG2_TF_neg/"
dir.create(myDir1)
mytemplate = "/lorax/chenglab/yanding/ATAC_seq_integration/Scripts/Template/HepG2/Template_for_counting_HepG2_TF_sig_neg.R"
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

}
length(bins)
[1] 895736

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
[1] 895736

rm(list=ls())
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/TF_counts_bins_logFC/"
myinf1 = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/TF_file_annotation_bigWig.txt"

files = list.files(mydir)
name = gsub(".Rda","",files)
name = gsub("_ENCFF172DEA_positive","",name)
name = gsub("_ENCFF172DEA_negative","",name)
name = unique(name)

tar = c("ENCFF285JUZ","ENCFF942IMP")
tag = which(name %in% tar)
name = name[-tag]

info = read.table(myinf1,sep="\t",quote=NULL,stringsAsFactors=F)
row.names(info) = info$file_accession
info = info[name,]
info = info[,c("target","file_accession")]

library(randomForest)
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Multi_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/RF_Uni_variable.R")
source("/lorax/chenglab/yanding/BRCA_PCR_RD/Scripts/Function/moveme.R")

AUC_score = list()

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_positive.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_negative.Rda")


tag_pos = sample(seq(1,895736),5000)
tag_neg = sample(seq(1,895736),5000)

rm(positive)
rm(negative)

res = matrix(0, 10000, nrow(info))
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
tag = c(rep(1,5000),rep(0,5000))
data = cbind(tag,res)
data = as.data.frame(data)
data$tag = factor(data$tag, levels =c(0,1))
fit = ROC(data)

1-fit[[3]]
[1] 0.8772474 0.8775325 0.8769620 0.8775456 0.8768321 0.8778768 0.8777907
 [8] 0.8767283 0.8765690 0.8774515
myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_TF_integration_AUC.Rda"
save.image(file = myoutf)

