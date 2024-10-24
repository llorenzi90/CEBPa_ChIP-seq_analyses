## ---------------------------
##
##
## Purpose of script:DiffBind with ER samples
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-05-09
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes: the problem is that I don't have good replicates 
##   
##
## ---------------------------

options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
## ---------------------------
library(DiffBind)
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/analyses/DiffBind/")
path_files=read.table("ER_bam_files.txt")

path_files$sample=gsub(".sorted.markedDups.proper_pairs.minq2.bam",".apr22",gsub(".sorted.markedDups.bam",".nov21",basename(path_files$V1)))
example_sample_sheet=read.csv("../sample_sheet_csuc-narrow.csv")
path_files$Condition=ifelse(grepl("p30",path_files$sample),"p30","p42")
path_files$Treatment=ifelse(grepl("LPS",path_files$sample),"LPS","UT")
path_files$Batch=ifelse(grepl("apr22",path_files$sample),"Apr22","Nov21")
path_files$Replicate=gsub("R","",ifelse(grepl("R3",path_files$sample),"R3",ifelse(grepl("R2",path_files$sample),"R2","R1")))

path_files$ID=paste(path_files$Condition,
                    path_files$Treatment,
                    path_files$Batch,sep = ".")

input_files=path_files[grepl("INP",path_files$sample),]
sample_files=path_files[!grepl("INP",path_files$sample),]
length(unique(input_files$ID))
length(unique(sample_files$ID))
table(input_files$ID%in%sample_files$ID)
table(sample_files$ID%in%input_files$ID)

sample_files$bamControl=input_files$V1[match(sample_files$ID,
                                             input_files$ID)]
sample_files$ControlID=input_files$sample[match(sample_files$ID,
                                             input_files$ID)]
sample_files$SampleID=paste(sample_files$Condition,
                            sample_files$Treatment,
                            sample_files$Replicate,
                            sample_files$Batch,sep = "_")

sample_files$Factor="ER"

sample_files$Peaks=paste0("/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/ER_merged_peaks/",sample_files$Condition,
                          "_",sample_files$Factor,"_",sample_files$Treatment,
                          "_merged_reps.bed")

sample_files$PeakCaller="bed"

colnames(sample_files)[1]="bamReads"
colnames(sample_files)
write.csv(sample_files[,c(10,3,4,11,5,6,1,8,9,12,13)],"sample_sheet_ER.csv",row.names = F)


##add inputs as samples to test how they cluster compared to bad ER samples

sf=sample_files[,c(10,3,4,11,5,6,1,8,9,12,13)]
colnames(sf)
input_files$SampleID=paste(input_files$Condition,
                            input_files$Treatment,
                            input_files$Batch,
                            "INP",sep="_")
input_files$Factor="ER"
colnames(input_files)[1]="bamReads"
input_files$bamControl=input_files$bamReads
input_files$ControlID=input_files$sample
input_files$Peaks=sf$Peaks[match(input_files$ID,
                                 sample_files$ID)]
input_files$PeakCaller="bed"

input_files$Condition=paste0(input_files$Condition,"INP")
sample_sheet_with_input=rbind(sf,input_files[,match(colnames(sf),colnames(input_files))])
write.csv(sample_sheet_with_input,"sample_sheet_wiht_input_ER.csv",row.names = F)

#generate also a sample sheet including K27 samples
#sample_sheet_K27

