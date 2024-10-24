## ---------------------------
##
##
## Purpose of script: generate input files for homer CEBPa ChiP-seq
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-03-04
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
## ---------------------------
datadir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/DESeq2_results/"
peaks_file <- "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/macs2_pooled_reps/merged_peaks/merged_peaks_EV.p30.p42_LPS.UT.sorted.bed.withcolnames.annotated.csv"
peaks <- read.csv(peaks_file)
setwd(datadir)
list.files()
outdir="bed_files_homer_input/"
pat="_DESeq2_res.noNAs.csv"
res_files <- list.files(pattern = pat)

# BED files should have at minimum 6 columns (separated by TABs, additional columns will be ignored)
# Column1: chromosome
# Column2: starting position
# Column3: ending position
# Column4: Unique Peak ID
# Column5: not used
# Column6: Strand (+/- or 0/1, where 0="+", 1="-")

colnames(peaks)

cols_to_keep <-c(1:3,19,8,5) 
head(peaks[,cols_to_keep])
peaks$strand <- "."
padjco <- 0.05
for (rf in res_files) {
  res <- read.csv(rf)
  comp <- gsub(pat,"",rf)
  DiffPeaks <- res[res$padj<=padjco,]
  
  upPeakIDs <- DiffPeaks$X[DiffPeaks$log2FoldChange>0]
  downPeakIDs <- DiffPeaks$X[DiffPeaks$log2FoldChange<0]
  backgroundPeakIDs <- res$X[res$padj>padjco]
  
  list_peak.sets <- list(upPeaks=upPeakIDs,
                         downPeaks=downPeakIDs,
                         backgroundPeaks=backgroundPeakIDs)
  
  for (peakset in names(list_peak.sets)) {
    peakids <- list_peak.sets[[peakset]]
    
    write.table(peaks[peaks$peakID%in%peakids,cols_to_keep],
                paste0(outdir,comp,"_",peakset,".bed"),
                col.names = F,row.names = F,quote = F, sep = "\t")
  }
 
}
