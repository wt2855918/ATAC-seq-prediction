########################################################
#GM12878
########################################################
[1]Histone download
rm(list=ls())
library(ENCODExplorer)
res <- queryEncode(organism = "Homo sapiens", assay = "ChIP-seq",
                      biosample_name = "GM12878", file_format = "bam",
                      fixed = TRUE)                     
tag = res$output_type  == "alignments"
res = res[tag,]
tag = grep("GRCh38",res$assembly)
res = res[tag,]
tag = grep("histone",res$investigated_as)
res = res[tag,]
#res = res[1,]
downloadEncode(res,dir ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/Histone_ChIP_seq/")

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/Histone_file_annotation.txt"
write.table(res,myoutf,sep="\t",quote=F)

[2]TF download
rm(list=ls())
library(ENCODExplorer)
res <- queryEncode(organism = "Homo sapiens", assay = "ChIP-seq",
                      biosample_name = "GM12878", file_format = "bam",
                      fixed = TRUE)                     
tag = res$output_type  == "alignments"
res = res[tag,]
tag = grep("GRCh38",res$assembly)
res = res[tag,]
tag = grep("transcription factor",res$investigated_as)
res = res[tag,]
#res = res[1,]
downloadEncode(res,dir ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq/")

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/TF_file_annotation.txt"
write.table(res,myoutf,sep="\t",quote=F)

#begin to examine the left over files
mydir = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq/"
files = list.files(mydir)
my_nam = gsub(".bam","",files)

tag = which(res$file_accession %in% my_nam)
res = res[-tag,]
downloadEncode(res,dir ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq/")

[3]Histone download by bigwig
rm(list=ls())
library(ENCODExplorer)
res <- queryEncode(organism = "Homo sapiens", assay = "ChIP-seq",
                      biosample_name = "GM12878", file_format = "bigWig",fuzzy = TRUE)                 


tag = res$output_type  == "fold change over control"
res = res[tag,]
tag = grep("hg19",res$assembly)
res = res[tag,]
tag = grep("histone",res$investigated_as)
res = res[tag,]
#res = res[1,]
downloadEncode(res,dir ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/Histone_ChIP_seq_bigWig/")

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/Histone_file_annotation_bigwig.txt"
write.table(res,myoutf,sep="\t",quote=F)

[4]TF download by bigwig
rm(list=ls())
library(ENCODExplorer)
res <- queryEncode(organism = "Homo sapiens", assay = "ChIP-seq",
                      biosample_name = "GM12878", file_format = "bigWig",
                      fixed = TRUE)                     
tag = res$output_type  == "fold change over control"
res = res[tag,]
tag = grep("hg19",res$assembly)
res = res[tag,]
tag = grep("transcription factor",res$investigated_as)
res = res[tag,]
#res = res[1,]
downloadEncode(res,dir ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/TF_ChIP_seq_bigWig/")

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/TF_file_annotation_bigWig.txt"
write.table(res,myoutf,sep="\t",quote=F)

[5]HIC data download
rm(list=ls())
library(ENCODExplorer)
res <- queryEncode(organism = "Homo sapiens", assay = "Hi-C",
                      biosample_name = "GM12878", file_format = "hic",
                      fixed = TRUE)                     

downloadEncode(res,dir ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/Hi_C_hic/")

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/GM12878/Hi_C_file_annotation_hic.txt"
write.table(res,myoutf,sep="\t",quote=F)

#batch download
/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/GM12878/



########################################################
#HepG2
########################################################
[1]Downd load TF hepG2
rm(list=ls())
library(ENCODExplorer)
res <- queryEncode(organism = "Homo sapiens", assay = "ChIP-seq",
                      biosample_name = "HepG2", file_format = "bigWig",
                      fixed = TRUE)                     
tag = res$output_type  == "fold change over control"
res = res[tag,]
tag = grep("hg19",res$assembly)
res = res[tag,]
tag = grep("transcription factor",res$investigated_as)
res = res[tag,]
#res = res[1,]
downloadEncode(res,dir ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/TF_ChIP_seq_bigWig/")

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/TF_file_annotation_bigWig.txt"
write.table(res,myoutf,sep="\t",quote=F)

[2]Downd load histone hepG2
rm(list=ls())
library(ENCODExplorer)
res <- queryEncode(organism = "Homo sapiens", assay = "ChIP-seq",
                      biosample_name = "HepG2", file_format = "bigWig",
                      fixed = TRUE)                     
tag = res$output_type  == "fold change over control"
res = res[tag,]
tag = grep("hg19",res$assembly)
res = res[tag,]
tag = grep("histone",res$investigated_as)
res = res[tag,]
#res = res[1,]
downloadEncode(res,dir ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/Histone_ChIP_seq_bigWig/")

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/Histone_file_annotation_bigWig.txt"
write.table(res,myoutf,sep="\t",quote=F)

[3]ECLIP data pos strand
rm(list=ls())
library(ENCODExplorer)
res <- queryEncode(organism = "Homo sapiens",  assay = "eCLIP",      
                      biosample_name = "HepG2", file_format = "bigWig",
                      fixed = TRUE)        
                      
tag = grep("hg19",res$assembly)
res = res[tag,]
tag = res$output_type == "plus strand signal of unique reads"
res = res[tag,]

dir.create("/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/RBP_eCLIP_seq_bigWig_pos/")
downloadEncode(res,dir ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/RBP_eCLIP_seq_bigWig_pos/")

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/RBP_file_annotation_bigWig_pos.txt"
write.table(res,myoutf,sep="\t",quote=F)

[4]ECLIP data neg strand
rm(list=ls())
library(ENCODExplorer)
res <- queryEncode(organism = "Homo sapiens",  assay = "eCLIP",      
                      biosample_name = "HepG2", file_format = "bigWig",
                      fixed = TRUE)        
                      
tag = grep("hg19",res$assembly)
res = res[tag,]
tag = res$output_type == "minus strand signal of unique reads"
res = res[tag,]

dir.create("/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/RBP_eCLIP_seq_bigWig_neg/")
downloadEncode(res,dir ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/HepG2/RBP_eCLIP_seq_bigWig_neg/")

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/HepG2/RBP_file_annotation_bigWig_neg.txt"
write.table(res,myoutf,sep="\t",quote=F)

########################################################
#mESC
########################################################
[1]Downd load HM mESC
rm(list=ls())
library(ENCODExplorer)
res <- queryEncode(organism = "Mus musculus", assay = "ChIP-seq",
                      biosample_name = "E14TG2a.4", file_format = "bigWig",
                      fixed = TRUE)                     
tag = res$output_type  == "fold change over control"
res = res[tag,]
tag = grep("mm10",res$assembly)
res = res[tag,]
tag = grep("histone",res$investigated_as)
res = res[tag,]
#res = res[1,]
downloadEncode(res,dir ="/lorax/chenglab/yanding/ATAC_seq_integration/Data/Raw_Data/E14TG2a/Histone_ChIP_seq_bigWig/")

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/Annotation/E14TG2a/Histone_file_annotation_bigWig.txt"
write.table(res,myoutf,sep="\t",quote=F)

info <- get_encode_df_full()
tag <- info$biosample_name == "GM12878"
info = info[tag,]

tag = grep("hg19",info$assembly)
info = info[tag,]
tag = info$file_format == "bigWig"
info = info[tag,]
tag = grep("histone",info$investigated_as)
info = info[tag,]