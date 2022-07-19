options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
## ---------------------------

library(bedtoolsr)
library(ggvenn)
source("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/cbPalette.R")
source("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/compute_plot_overlap.R")  
#Overlap between replicates ER UT
#new replicate
p30ER.UT.R1=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/ER/peaks/apr22/p30_UT_R1_ER.sorted.markedDups.proper_pairs.minq2.bam.NoModel_q_0.05_peaks.narrowPeak")
#old replicate
p30ER.UT.old=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/ER/peaks/nov21/p30_ER.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")

#compute merging/overlap
ol=compute_plot_overlap(list(oldRep_nov21=p30ER.UT.old,newRep_apr22=p30ER.UT.R1),tit="p30 ER overlap replicates")

#check scores for merged peaks
View(ol[[1]])

library(ggplot2)
dtp=ol[[1]]

ggplot(dtp,aes(x=V6,y=10^-(V7))) + geom_boxplot() + theme_classic()

pdf("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/peaksets_overlaps/overlap_old_new_p30_persetScores.pdf")
ggplot(dtp,aes(x=V6,y=log10(V7))) + geom_boxplot() + theme_classic()+
  ylab("log10(score)") + xlab("")

dev.off()
tapply(dtp$V7, dtp$V6, summary)

#apply different cutoffs
list_of_files=list(oldRep_nov21=p30ER.UT.old,newRep_apr22=p30ER.UT.R1)
for (co in c(1:20)) {
  filt_list=lapply(list_of_files, function(x)return(x[x$V9>=co,]))
  print(lapply(filt_list, nrow))
  o=compute_plot_overlap(filt_list)
  print(table(o[[1]]$V6))
}

plot(density(log2(dtp$V7)))
dtp$width=dtp$V3 - dtp$V2

pdf("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/peaksets_overlaps/overlap_old_new_p30_persetWidths.pdf")

ggplot(dtp,aes(x=V6,y=log10(width))) + geom_boxplot() + theme_classic()+
  ylab("log10(width)") + xlab("")
dev.off()

tapply(dtp$width, dtp$V6, summary)
