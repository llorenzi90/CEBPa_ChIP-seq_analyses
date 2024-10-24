setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/comparison_Jakobsen_et_al_2019_peaks/")
peaks_table=readxl::read_xlsx("aaw4304_data_file_s1.xlsx",sheet = 1)
table(peaks_table$class)
peaks_table$name=paste0(peaks_table$chr,":",peaks_table$start,"-",peaks_table$end)
peaks_table$score=as.numeric(factor(peaks_table$class,levels=c("p30","p42","common")))
peaks_table$score=peaks_table$class
peaks_table$strand="."
#bed format: chr,start,end,name,score,strand
write.table(peaks_table[,c(1:3,6:8)],"Jakobsen2019_CEBPa_peaks.bed",quote = F,col.names = F,row.names = F,sep = "\t")

