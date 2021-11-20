[1]sequence define positive sequence
rm(list=ls())
library("Biostrings")
library("reticulate")
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1/ENCFF172DEA_positive_peak.fa"
fastaFile <- readDNAStringSet(myinf)

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_positive.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_tag.Rda")

#restrain the shared positive
ID = paste0(positive$V1,"_",positive$V2,"_",positive$V3)
row.names(positive) = ID
com = intersect(ID, names(fastaFile))
positive = positive[com,]
fastaFile = fastaFile[com]
ID = com

#pick up the label
total_seq = seq(1,length(com),1)
training_pos_tag = sample(seq(1,length(com)),20000)
validation_pos_pool = setdiff(total_seq, training_pos_tag)
validation_pos_tag = sample(validation_pos_pool,5000)

#begin to get the sequence file

ID_for_training = ID[training_pos_tag]

training_fasta = fastaFile[ID_for_training]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1_training_and_validation/ENCFF172DEA_positive_training.fa"
writeXStringSet(training_fasta, myoutf)

ID_for_validation = ID[validation_pos_tag]
validation_fasta = fastaFile[ID_for_validation]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1_training_and_validation/ENCFF172DEA_positive_validation.fa"
writeXStringSet(validation_fasta, myoutf)

#begin to get the bed files

training_pos = positive[training_pos_tag,]
validation_pos = positive[validation_pos_tag,]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/training_pos_bed_file.txt"
write.table(training_pos, myoutf, col.names=F,row.names=F,quote=F)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/validation_pos_bed_file.txt"
write.table(validation_pos, myoutf, col.names=F,row.names=F,quote=F)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_pos.Rda"
save.image(file = myoutf)

[3]sequence define negative sequence
rm(list=ls())
library("Biostrings")
library("reticulate")
myinf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1/ENCFF172DEA_negative_peak.fa"
fastaFile <- readDNAStringSet(myinf)

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_peak_processed_bed/Batch1/ENCFF172DEA_negative.Rda")
load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_tag.Rda")

#restrain the peak that can be plot the sequence
ID = paste0(negative$V1,"_",negative$V2,"_",negative$V3)
row.names(negative) = ID
com = intersect(ID, names(fastaFile))
negative = negative[com,]
fastaFile = fastaFile[com]
ID = com

#pick up the label
total_seq = seq(1,length(com),1)
training_neg_tag = sample(seq(1,length(com)),20000)
validation_neg_pool = setdiff(total_seq, training_neg_tag)
validation_neg_tag = sample(validation_neg_pool,5000)

#begin to get sequence file
ID_for_training = ID[training_neg_tag]

training_fasta = fastaFile[ID_for_training]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1_training_and_validation/ENCFF172DEA_negative_training.fa"
writeXStringSet(training_fasta, myoutf)

ID_for_validation = ID[validation_neg_tag]
validation_fasta = fastaFile[ID_for_validation]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/ATAC_seq_bin_sequence/Batch_1_training_and_validation/ENCFF172DEA_negative_validation.fa"
writeXStringSet(validation_fasta, myoutf)

#begin to get the bed files
training_neg = negative[training_neg_tag,]
validation_neg = negative[validation_neg_tag,]

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/training_neg_bed_file.txt"
write.table(training_neg, myoutf, col.names=F,row.names=F,quote=F)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/validation_neg_bed_file.txt"
write.table(validation_neg, myoutf, col.names=F,row.names=F,quote=F)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/GM12878_training_and_validation_neg.Rda"
save.image(file = myoutf)

[4]bed profile prepration

cd /lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/

cut -f 1,2,3 training_pos_bed_file.txt > training_pos.bed
cut -f 1,2,3 validation_pos_bed_file.txt > validation_pos.bed

cut -f 1,2,3 training_neg_bed_file.txt > training_neg.bed
cut -f 1,2,3 validation_neg_bed_file.txt > validation_neg.bed

