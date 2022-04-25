## ---------------------------
##
##
## Purpose of script: add annotation to DESeq2 results
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-03-29
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
## ---------------------------options(scipen = 999) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)
# library(ChIPseeker)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## ---------------------------
peakfile="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/macs2_pooled_reps/merged_peaks/merged_peaks_EV.p30.p42_LPS.UT.sorted.bed.withcolnames"
peakAnnodf <- read.csv(paste0(peakfile,".annotated.csv"))
#read in DESeq2 results
deseq2resdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/DESeq2_results/"
setwd(deseq2resdir)
files=list.files(pattern = "_DESeq2_res.noNAs.csv")
annot_files=list.files("annotated_results/")
paste0(gsub(".csv","",files),".annotated.csv")%in%annot_files
files=files[!paste0(gsub(".csv","",files),".annotated.csv")%in%annot_files]

for (fi in files) {
  deseq2res=read.csv(fi)
  #add annotation
  deseq2res <- cbind(deseq2res,peakAnnodf[match(deseq2res$X,peakAnnodf$peakID),c(6:18)])
  #write annotated results
  write.csv(deseq2res, paste0("annotated_results/",gsub(".csv","",fi),".annotated.csv"),row.names = F)
}

