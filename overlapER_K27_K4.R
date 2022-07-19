## ---------------------------
##
##
## Purpose of script:compare ER peaks with k27 and k4 peaks
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-05-10
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
source("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/compute_plot_overlap.R")
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/")
#ER peaks
p30ER=read.table("data/ER/peaks/nov21/p30_ER.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")
p42ER=read.table("data/ER/peaks/nov21/p42_ER.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")


#K27 peaks 
#p30:
p30K27_list=lapply(1:3, function(x)read.table(paste0("data/K27/peaks/p30_UT_R",x,".sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")))
names(p30K27_list)=paste0("p30.K27.R",1:3)
p30K27_list$p30.K27.R1

compute_plot_overlap(p30K27_list,tit = "p30-K27")
compute_plot_overlap(c(p30K27_list,list(p30.ER=p30ER)))

#K27 pooled reps
p30K27_pooled=read.table("data/K27/peaks/pooled_reps/p30_UT.pooled.NoModel_q_0.05_peaks.narrowPeak")

compute_plot_overlap(list(ER=p30ER,K27_pooled=p30K27_pooled),tit="p30 ER vs K27 pooled reps")

#K27 pooled reps R1 and R3
p30K27_pooled_R1R3=read.table("data/K27/peaks/pooled_reps/p30_UT_R1_R3_pooled_reps.pooled.NoModel_q_0.05_peaks.narrowPeak")

compute_plot_overlap(list(ER=p30ER,K27_pooled.R1R3=p30K27pooled_R1R3),tit="p30 ER vs K27 pooled reps")


#p42
p42K27_list=lapply(1:2, function(x)read.table(paste0("data/K27/peaks/p42_UT_R",x,".sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")))
names(p42K27_list)=paste0("p42.K27.R",1:2)
p42K27_list$p42.K27.R1

compute_plot_overlap(p42K27_list,tit = "p42-K27")
compute_plot_overlap(c(p42K27_list,list(p42.ER=p42ER)))

#K27 pooled reps
p42K27_pooled=read.table("data/K27/peaks/pooled_reps/p42_UT.pooled.NoModel_q_0.05_peaks.narrowPeak")

compute_plot_overlap(list(ER=p42ER,K27_pooled=p42K27_pooled),tit="p42 ER vs K27 pooled reps")


#the same for K4...
p30K4_list=lapply(1:2, function(x)read.table(paste0("data/K4/peaks/p30_UT_R",x,"_K4.sorted.markedDups.proper_pairs.minq2.bam.NoModel_q_0.05_peaks.narrowPeak")))
names(p30K4_list)=paste0("p30.K4.R",1:2)

p30K4_pooled=read.table("data/K4/peaks/pooled_reps/p30_UT_K4.pooled.NoModel_q_0.05_peaks.narrowPeak")


p42K4_list=lapply(1:2, function(x)read.table(paste0("data/K4/peaks/p42_UT_R",x,"_K4.sorted.markedDups.proper_pairs.minq2.bam.NoModel_q_0.05_peaks.narrowPeak")))
names(p42K4_list)=paste0("p42.K4.R",1:2)

p42K4_pooled=read.table("data/K4/peaks/pooled_reps/p42_UT_K4.pooled.NoModel_q_0.05_peaks.narrowPeak")


#######Part A#############
#Do it automatically: for each vector (p30 and p42) ER comparison do:
# compare against each factor (k27 and K4):
# check how many replicates for UT
# load individual replicates peaks into a list
#1) compute and save plot overlap with all replicates 
#2) compute and save plot overlap with merged replicates
#3) compute and save plot overlap with pooled replicates

vectors=c("p30","p42")

#read all peaks files that we will use
all_peaks=c(list(p30.ER=p30ER,p42.ER=p42ER),list(p30.K27.pooled=p30K27_pooled,
                                                 p42.K27.pooled=p42K27_pooled),
            p30K27_list,p42K27_list,list(p30.K4.pooled=p30K4_pooled,
                                         p42.K4.pooled=p42K4_pooled),
            p30K4_list,p42K4_list)

setwd("analyses/ER/")
dir.create("Overlap_K27_K4")
setwd("Overlap_K27_K4/")

for (fac in c("K27","K4")) {
  for (vec in vectors) {
    ER=all_peaks[paste(vec,"ER",sep = ".")]
    li=all_peaks[grep(paste0(vec,".",fac,".R"),names(all_peaks))]
    pool=all_peaks[paste(vec,fac,"pooled",sep = ".")]
    
    #overlap with all replicates
    c=compute_plot_overlap(c(ER,li),tit = paste0(vec," ER vs ",fac," individual reps"))

    pdf(paste0(vec,"_ER_vs_",fac,"_individual_reps.pdf"))
    print(c[[2]])
    dev.off()
    
    #overlap with merged replicates
    mergedreps=list(compute_plot_overlap(li)[[1]])
    names(mergedreps)=paste0(vec,".",fac,".mergedreps")
    c = compute_plot_overlap(c(ER,mergedreps),tit = paste0(vec," ER vs ",fac," merged reps"))

    pdf(paste0(vec,"_ER_vs_",fac,"_merged_reps.pdf"))
    print(c[[2]])
    dev.off()
    
    #overlap with pooled replicates
    c=compute_plot_overlap(c(ER,pool),tit = paste0(vec," ER vs ",fac," pooled reps"))
    
    pdf(paste0(vec,"_ER_vs_",fac,"_pooled_reps.pdf"))
    print(c[[2]])
    dev.off()
  }
}


