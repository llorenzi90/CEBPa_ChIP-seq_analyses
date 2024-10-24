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
## Notes: First version of this script was WRONG!!! I used the human annotation instead of mouse annotation!!
##   This has been corrected now (20 may 22)
##
## ---------------------------

options(scipen = 999) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
## ---------------------------

peakfile <- "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/macs2_pooled_reps/merged_peaks/merged_peaks_EV.p30.p42_LPS.UT.sorted.bed"
peakfile <- "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/K27/peaks/merged_peaks/merged_peaks_EV.p30.p42_LPS.UT.sorted.bed"

txdb <- TxDb.Mmusculus.UCSC.mm39.refGene
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
#gene_translation_file <- "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/biomart_tables/Human_Ensembl_genes_104_GRCh38.p13_to_NCBI.tsv"
#that was wrong!! human
library(data.table)
#gene_translation <- fread(gene_translation_file)
# table(peakAnnodf$geneId%in%gene_translation$`NCBI gene (formerly Entrezgene) ID`)
# peakAnnodf$gene_name <- gene_translation$`Gene name`[match(peakAnnodf$geneId,
#                                                            gene_translation$`NCBI gene (formerly Entrezgene) ID`)]
# peakAnnodf$peakID=paste0(peakAnnodf$seqnames,":",peakAnnodf$start - 1,"-",peakAnnodf$end)

gene_translation_file <- "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/mouse/conversion_table_ENSEMBL_NCBI_ids.csv"
gene_translation <- fread(gene_translation_file)

peakAnnodf <- as.data.frame(peakAnno)
table(peakAnnodf$geneId%in%gene_translation$ncbi_GeneID)
peakAnnodf$gene_name <- gene_translation$gene_name[match(peakAnnodf$geneId,
                                                         gene_translation$ncbi_GeneID)]
peakAnnodf$peakID=paste0(peakAnnodf$seqnames,":",peakAnnodf$start -1 ,"-",peakAnnodf$end)
write.csv(peakAnnodf,paste0(peakfile,".annotated.csv"),row.names = F)
#read in DESeq2 results
deseq2resdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/K27/macs2_merged_peaks_counts/DESeq2_results/"
setwd(deseq2resdir)
dir.create("annotated_results")
files=list.files(pattern = "_DESeq2_res.noNAs.csv")
for (fi in files) {
  deseq2res=read.csv(fi)
  #add annotation
  deseq2res <- cbind(deseq2res,peakAnnodf[match(deseq2res$X,peakAnnodf$peakID),c(6:18)])
  #write annotated results
  write.csv(deseq2res, paste0("annotated_results/",gsub(".csv","",fi),".annotated.csv"),row.names = F)
}

