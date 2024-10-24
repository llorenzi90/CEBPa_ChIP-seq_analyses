## ---------------------------
##
##
## Purpose of script: find overlap between ChIP-seq peaks
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-05-05
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

library(bedtoolsr)
library(ggvenn)
source("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/cbPalette.R")
compute_plot_overlap <- function(mylist,tit=""){
  library(bedtoolsr)
  library(ggvenn)
  
  le=length(mylist)
  if(is.null(names(mylist))){
    names(mylist)=paste0("file",seq(1:le))
  }
  na=names(mylist)
  
  #check columns are the same (10 columns):
  mylist=lapply(mylist,function(x){
    if(ncol(x)!=10){
      y=x[,c(1:ncol(x),rep(ncol(x),10-ncol(x)))]
      colnames(y)=paste0("V",seq(1:10))
      return(y)
    }else{colnames(x)=paste0("V",seq(1:10))
    return(x)}
  })
  mylist <- lapply(names(mylist),function(x){
    y=mylist[[x]]
    y$V4=x
    return(y)
  })
  ca=do.call("rbind",mylist)
  ca=ca[order(ca$V1,ca$V2),]
  merg=bt.merge(ca,c = "4,4,4,9",d=-1,o = "count,count_distinct,distinct,mean") 
  
  dtv=lapply(na, function(x)return(rownames(merg)[grepl(x,merg$V6)]))
  
  names(dtv)=na
  si=8
  if(tit==""){tit=paste(na,collapse = "_vs_")}
  g <- ggvenn(
    dtv, 
    fill_color = cbPalette[1:length(na)], #change this to set other colors 
    #fill_color = c( "#EFC000FF", "#868686FF"),
    stroke_size = 0.5, set_name_size = si,
    show_percentage = F,
    text_size = si
  ) + ggtitle(tit) + theme(plot.title =  element_text(size=si*2.5))
  
  
  
  return(list(merg,g))
}


setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/Venn_diagrams/peaksets_overlaps/comparisons_ER/")

#Overlap between replicates ER LPS p30
p30ER.LPS.R1=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/may22/p30_LPS_R1_ER.sorted.markedDups.proper_pairs.minq2.bam.NoModel_q_0.05_peaks.narrowPeak")
p30ER.LPS.R2=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/may22/p30_LPS_R2_ER.sorted.markedDups.proper_pairs.minq2.bam.NoModel_q_0.05_peaks.narrowPeak")
p30ER.LPS.R3=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/may22/p30_LPS_R3_ER.sorted.markedDups.proper_pairs.minq2.bam.NoModel_q_0.05_peaks.narrowPeak")

p30_ER_LPS_reps=list(p30ER.LPS.R1=p30ER.LPS.R1,
                     p30ER.LPS.R2=p30ER.LPS.R2,
                     p30ER.LPS.R3=p30ER.LPS.R3)
pdf("p30_LPS_ER_overlap_between_replicates.pdf")
print(compute_plot_overlap(p30_ER_LPS_reps)[[2]])
dev.off()

#Overlap between pooled and merged

#1 extract merging of replicates from call of compute_plot_overlap
p30_LPS_ER_merged_reps=compute_plot_overlap(p30_ER_LPS_reps)[[1]]
#2 compute overlap between merged and pooled peaks
p30_LPS_ER.pooled=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/macs2_pooled_reps/may22/p30_LPS_ER.pooled.NoModel_q_0.05_peaks.narrowPeak")

pdf("p30_LPS_ER_merged_vs_pooled.pdf")
compute_plot_overlap(list(p30_LPS_ER_merged_reps=p30_LPS_ER_merged_reps,
                          p30_LPS_ER_pooled=p30_LPS_ER.pooled))
dev.off()

#Overlap between replicates ER UT
#new replicate
p30ER.UT.R1=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/may22/p30_UT_R1_ER.sorted.markedDups.proper_pairs.minq2.bam.NoModel_q_0.05_peaks.narrowPeak")
#old replicate
p30ER.UT.old=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/p30_ER.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")

pdf("p30_UT_ER_overlap_between_replicates_old_new.pdf")
compute_plot_overlap(list(oldRep_nov21=p30ER.UT.old,newRep_apr22=p30ER.UT.R1),tit="p30 ER overlap replicates")
dev.off()

#compute merging between new and old replicates
p30_ER_mergedreps=compute_plot_overlap(list(oldRep_nov21=p30ER.UT.old,newRep_apr22=p30ER.UT.R1))[[1]]

#Compute overlap between merged peaks both LPS vs ER
pdf("overlap_p30_ER_LPS_vs_UT_merged_reps.pdf")
compute_plot_overlap(list(p30_LPS_mergedreps=p30_LPS_ER_merged_reps,
                     p30_UT_mergedreps=p30_ER_mergedreps),tit="p30 ER LPS vs UT")
dev.off()

#Still to do: check p-values of each set for each comparison, to answer the question, 
#are those peaks that are only found in one condition (espeacially the old p30ERUT that has a lot of peaks)
#less confident than the ones that overlap? are between them very confident peaks


#Do the same for p42:
#...
#Overlap between replicates ER LPS p42
p42ER.LPS.R1=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/may22/p42_LPS_R1_ER.sorted.markedDups.proper_pairs.minq2.bam.NoModel_q_0.05_peaks.narrowPeak")
p42ER.LPS.R2=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/may22/p42_LPS_R2_ER.sorted.markedDups.proper_pairs.minq2.bam.NoModel_q_0.05_peaks.narrowPeak")


p42_ER_LPS_reps=list(p42ER.LPS.R1=p42ER.LPS.R1,
                     p42ER.LPS.R2=p42ER.LPS.R2)
pdf("p42_LPS_ER_overlap_between_replicates.pdf")
print(compute_plot_overlap(p42_ER_LPS_reps)[[2]])
dev.off()

#Overlap between pooled and merged

#1 extract merging of replicates from call of compute_plot_overlap
p42_LPS_ER_merged_reps=compute_plot_overlap(p42_ER_LPS_reps)[[1]]
#2 compute overlap between merged and pooled peaks
p42_LPS_ER.pooled=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/macs2_pooled_reps/may22/p42_LPS_ER.pooled.NoModel_q_0.05_peaks.narrowPeak")

pdf("p42_LPS_ER_merged_vs_pooled.pdf")
compute_plot_overlap(list(p42_LPS_ER_merged_reps=p42_LPS_ER_merged_reps,
                          p42_LPS_ER_pooled=p42_LPS_ER.pooled))
dev.off()


#Overlap between replicates ER UT
#new replicate
p42ER.UT.R1=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/may22/p42_UT_R2_ER.sorted.markedDups.proper_pairs.minq2.bam.NoModel_q_0.05_peaks.narrowPeak")
#old replicate
p42ER.UT.old=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/p42_ER.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")

pdf("p42_UT_ER_overlap_between_replicates_old_new.pdf")
compute_plot_overlap(list(oldRep_nov21=p42ER.UT.old,newRep_apr22=p42ER.UT.R1),tit="p42 ER overlap replicates")
dev.off()

#compute merging between new and old replicates
p42_ER_mergedreps=compute_plot_overlap(list(oldRep_nov21=p42ER.UT.old,newRep_apr22=p42ER.UT.R1))[[1]]

#Compute overlap between merged peaks both LPS vs ER
pdf("overlap_p42_ER_LPS_vs_UT_merged_reps.pdf")
compute_plot_overlap(list(p42_LPS_mergedreps=p42_LPS_ER_merged_reps,
                          p42_UT_mergedreps=p42_ER_mergedreps),tit="p42 ER LPS vs UT")
dev.off()



#p30 vs p42 
#nov21
pdf("ER_UT_p30_vs_p42_nov21.pdf")
compute_plot_overlap(list(p30UT=p30ER.UT.old,p42UT=p42ER.UT.old),tit="ER UT peaks p30 vs p42")
dev.off()

#apr22
pdf("ER_UT_p30_vs_p42_apr22.pdf")
compute_plot_overlap(list(p30UT=p30ER.UT.R1,p42UT=p42ER.UT.R1),tit="ER UT peaks p30 vs p42")
dev.off()

#merged reps old_new
#Compute overlap between merged peaks both LPS vs ER
pdf("overlap_ER_LPS_p30_vs_p42_merged_reps.pdf")
compute_plot_overlap(list(p30_LPS_mergedreps=p30_LPS_ER_merged_reps,
                          p42_LPS_mergedreps=p42_LPS_ER_merged_reps
                          ),tit="ER LPS p30 vs p42")
dev.off()

pdf("overlap_ER_UT_p30_vs_p42_merged_reps.pdf")
compute_plot_overlap(list(p30_UT_mergedreps=p30_ER_mergedreps,
                          p42_UT_mergedreps=p42_ER_mergedreps
),tit="ER UT p30 vs p42")
dev.off()

#overlap between all 4 UT ER

pdf("ER_UT_p30_vs_p42_individual_reps.pdf")
compute_plot_overlap(list(p30.old=p30ER.UT.old,
                          p30.new=p30ER.UT.R1,
                          p42.old=p42ER.UT.old,
                          p42.new=p42ER.UT.R1),tit="ER-UT")
dev.off()

#write merged peaks files
dir.create("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/data/merged_peaks")
write.table(p30_LPS_ER_merged_reps[,c(1:4,7,6)],"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/data/merged_peaks/p30_ER_LPS_merged_reps.bed",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(p30_ER_mergedreps[,c(1:4,7,6)],"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/data/merged_peaks/p30_ER_UT_merged_reps.bed",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(p42_LPS_ER_merged_reps[,c(1:4,7,6)],"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/data/merged_peaks/p42_ER_LPS_merged_reps.bed",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(p42_ER_mergedreps[,c(1:4,7,6)],"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/data/merged_peaks/p42_ER_UT_merged_reps.bed",row.names = F,col.names = F,quote = F,sep = "\t")


#overlap between K27 p30 UT

setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/")
p30_K27_UT_R1=read.table("p30_UT_R1.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")
p30_K27_UT_R2=read.table("p30_UT_R2.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")
p30_K27_UT_R3=read.table("p30_UT_R3.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")

overlap_p30_K27_UT=compute_plot_overlap(list(R1=p30_K27_UT_R1,
                                             R2=p30_K27_UT_R2,
                                             R3=p30_K27_UT_R3),tit="p30_K27_UT")
overlap_p30_K27_UT[[2]]


p30_K27_LPS_R1=read.table("p30_LPS_R1.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")
p30_K27_LPS_R2=read.table("p30_LPS_R2.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")
p30_K27_LPS_R3=read.table("p30_LPS_R3.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")

overlap_p30_K27_LPS=compute_plot_overlap(list(R1=p30_K27_LPS_R1,
                                             R2=p30_K27_LPS_R2,
                                             R3=p30_K27_LPS_R3),tit="p30_K27_LPS")
overlap_p30_K27_LPS[[2]]




p42_K27_UT_R1=read.table("p42_UT_R1.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")
p42_K27_UT_R2=read.table("p42_UT_R2.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")

overlap_p42_K27_UT=compute_plot_overlap(list(R1=p42_K27_UT_R1,
                                             R2=p42_K27_UT_R2),tit="p42_K27_UT")
overlap_p42_K27_UT[[2]]


p42_K27_LPS_R1=read.table("p42_LPS_R1.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")
p42_K27_LPS_R2=read.table("p42_LPS_R2.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")

overlap_p42_K27_LPS=compute_plot_overlap(list(R1=p42_K27_LPS_R1,
                                              R2=p42_K27_LPS_R2),tit="p42_K27_LPS")
overlap_p42_K27_LPS[[2]]



all_ER_merged=compute_plot_overlap(list(p30.LPS.ER=p30_LPS_ER_merged_reps,
                          p30.UT.ER=p30_ER_mergedreps,
                          p42.LPS.ER=p42_LPS_ER_merged_reps,
                          p42.UT.ER=p42_ER_mergedreps))
all_ER_merged= all_ER_merged[[1]]
write.table(all_ER_merged[,c(1:4,7,6)], "~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/data/merged_peaks/all_ER_mereged_peaks.bed", quote = F, row.names = F,col.names = F,
            sep = "\t")


