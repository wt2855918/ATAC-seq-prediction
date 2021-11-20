[1]Begin to annotate the differen data in GM12878 ATAC-seq
rm(list=ls())

#define the convert function
peakDF2GRanges <- function(peak.df) {
    peak.gr=GRanges(seqnames=peak.df[,1],
        ranges=IRanges(peak.df[,2], peak.df[,3]))
    cn <- colnames(peak.df)
    if (length(cn) > 3) {
        for (i in 4:length(cn)) {
            mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
        }
    }
    return(peak.gr)
}


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(GenomicRanges)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/Integration/HepG2_TF_integration_Core.Rda")
tag = res2$tag == "1"
res2 = res2[tag==1,]

chrom = as.vector(sapply(row.names(res2), function(x) strsplit(x,"_")[[1]][1]))
start = as.vector(sapply(row.names(res2), function(x) strsplit(x,"_")[[1]][2]))
end = as.vector(sapply(row.names(res2), function(x) strsplit(x,"_")[[1]][3]))

res = cbind(chrom,start,end)
res = as.data.frame(res)
res$start = as.numeric(as.vector(res$start))
res$end = as.numeric(as.vector(res$end))
res = peakDF2GRanges(res)
peakAnno <- annotatePeak(res, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_8/HepG2_ATAC_annotation.pdf"
pdf(myoutf,width=5,height=5)
plotAnnoPie(peakAnno)
dev.off()

info = peakAnno@anno
info = as.data.frame(info)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/positive_region_annotation.txt"
write.table(info,myoutf,sep="\t",quote=F)

[2]Begin to annotate the differen data in GM12878 DNase-seq
rm(list=ls())

#define the convert function
peakDF2GRanges <- function(peak.df) {
    peak.gr=GRanges(seqnames=peak.df[,1],
        ranges=IRanges(peak.df[,2], peak.df[,3]))
    cn <- colnames(peak.df)
    if (length(cn) > 3) {
        for (i in 4:length(cn)) {
            mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
        }
    }
    return(peak.gr)
}

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(GenomicRanges)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/DNase_seq_peak_processed_bed/ENCFF422EDI_positive.Rda")
res2 = positive
res2$V1 = as.vector(res2$V1)
res2$V2 = as.numeric(as.vector(res2$V2))
res2$V3 = as.numeric(as.vector(res2$V3))

res = peakDF2GRanges(res2)
peakAnno <- annotatePeak(res, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_8/HepG2_DNase_annotation.pdf"
pdf(myoutf,width=5,height=5)
plotAnnoPie(peakAnno)
dev.off()

info = peakAnno@anno
info = as.data.frame(info)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/positive_DNase_region_annotation.txt"
write.table(info,myoutf,sep="\t",quote=F)

[3]Beign to annotate the different data in GM12878 FAIR-seq
rm(list=ls())

#define the convert function
peakDF2GRanges <- function(peak.df) {
    peak.gr=GRanges(seqnames=peak.df[,1],
        ranges=IRanges(peak.df[,2], peak.df[,3]))
    cn <- colnames(peak.df)
    if (length(cn) > 3) {
        for (i in 4:length(cn)) {
            mcols(peak.gr)[[cn[i]]] <- peak.df[, cn[i]]
        }
    }
    return(peak.gr)
}
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(GenomicRanges)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

load("/lorax/chenglab/yanding/ATAC_seq_integration/Data/HepG2/FAIRE_seq_peak_processed_bed/ENCFF001UYN_positive.Rda")
res2 = positive

res2 = positive
res2$V1 = as.vector(res2$V1)
res2$V2 = as.numeric(as.vector(res2$V2))
res2$V3 = as.numeric(as.vector(res2$V3))

res = peakDF2GRanges(res2)
peakAnno <- annotatePeak(res, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Figures/MS_Figure/Figure_8/HepG2_FAIRE_annotation.pdf"
pdf(myoutf,width=5,height=5)
plotAnnoPie(peakAnno)
dev.off()

info = peakAnno@anno
info = as.data.frame(info)

myoutf = "/lorax/chenglab/yanding/ATAC_seq_integration/Data/GM12878/positive_DNase_region_annotation.txt"
write.table(info,myoutf,sep="\t",quote=F)
