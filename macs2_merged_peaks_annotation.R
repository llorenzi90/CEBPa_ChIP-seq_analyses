## ---------------------------
##
##
## Purpose of script:peak annotation CEBPa ChIP-peaks from merged macs2 peaks across samples (from pooled reps) 
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-03-02
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

options(scipen = 999) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## ---------------------------

peakfile <- "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/macs2_pooled_reps/merged_peaks/merged_peaks_EV.p30.p42_LPS.UT.sorted.bed"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks=read.table(peakfile)
colnames(peaks) <- c( "seqnames",
                      "start",
                      
                      "end",
                      "nsamples",
                      "samples",
                      "meanscore")
write.table(peaks,paste0(peakfile,".withcolnames"),row.names = F,quote = F,sep = "\t")
peakfile <- paste0(peakfile,".withcolnames")

peakAnno <- annotatePeak(peakfile, tssRegion=c(-3000, 3000), TxDb=txdb)
View(as.data.frame(peakAnno))

#add some more info: translate geneID, add p-values
gene_translation_file <- "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/biomart_tables/Human_Ensembl_genes_104_GRCh38.p13_to_NCBI.tsv"
library(data.table)
gene_translation <- fread(gene_translation_file)
peakAnnodf <- as.data.frame(peakAnno)
table(peakAnnodf$geneId%in%gene_translation$`NCBI gene (formerly Entrezgene) ID`)
peakAnnodf$gene_name <- gene_translation$`Gene name`[match(peakAnnodf$geneId,
                                                           gene_translation$`NCBI gene (formerly Entrezgene) ID`)]
peakAnnodf$peakID=paste0(peakAnnodf$seqnames,":",peakAnnodf$start - 1,"-",peakAnnodf$end)
write.csv(peakAnnodf,paste0(peakfile,".annotated.csv"),row.names = F)
#read in DESeq2 results
deseq2resdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/DESeq2_results/"
setwd(deseq2resdir)
files=list.files(pattern = "_DESeq2_res.noNAs.csv")
for (fi in files) {
  deseq2res=read.csv(fi)
  #add annotation
  deseq2res <- cbind(deseq2res,peakAnnodf[match(deseq2res$X,peakAnnodf$peakID),c(6:18)])
  #write annotated results
  write.csv(deseq2res, paste0("annotated_results/",gsub(".csv","",fi),".annotated.csv"),row.names = F)
}

