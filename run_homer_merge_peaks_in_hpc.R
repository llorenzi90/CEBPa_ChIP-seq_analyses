#notes:
#homer format:
# 1. Merged Peak name (will start with "Merged-")
# 2. chromosome
# 3. start (average from merged peaks)
# 4. end (average from merged peaks)
# 5. strand
# 6. Average peak score (actually, the average of the original values in column 6 of the peak files - or column 5 of BED files)
# 7. Original peak files contributing to the merged peak
# 8. Total number of peaks merged (occasionally more than one peak from a single file will be merged if the peaks are within the specify distance or two or more peaks from one file overlap with the same single peak(s) from another file)



#####Run on cluster ####
#before running type in login hpc:

#module load conda
#conda activate csaw (homer is installed in this environment)
#R
#> #(type the following commands:)

setwd("/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples")
getwd()
rps=read.table("reps_per_sample.txt")
rps
# system("module load conda/current")
# system("conda activate csaw")
setwd("merged_peaks")
for (i in 1:nrow(rps)) {
  samp=gsub("\\./","",rps$V2[i])
  nrep=rps$V1[i]
  if(nrep==1){#if no reps convert macs2 peak file to homer format
    cat("sample:")
    system(paste0("ls .",rps$V2[i],"*/*narrowPeak"))
    tmp=read.table(system(paste0("ls ../",samp,"*/*narrowPeak"),intern = T))
    write.table(data.frame(tmp$V4,tmp$V1,tmp$V2+1,tmp$V3,"+",tmp$V5,gsub("_peak_[0-9]*","_peaks.narrowPeak",tmp$V4),"1"),paste0(samp,".homer.peaks"),col.names = F,row.names = F,quote = F,sep = "\t")
  } else{ cat("samples:")
  system(paste0("ls .",rps$V2[i],"*/*narrowPeak"))
  system(paste0("mergePeaks -venn ",samp,".venn -matrix ",samp,".matrix -d given .",rps$V2[i],"*/*narrowPeak -prefix ", samp))
  system(paste0("mv ",samp,paste(paste0("*R",seq_len(nrep)),collapse = ""),"*narrowPeak ",samp,".mergedReps.homer.peaks"))
  }
}


#all pairwise comparisons:
compars_mat=combn(gsub("\\./","",rps$V2),2)
for (i in 1:ncol(compars_mat)) {
  samps=compars_mat[,i]
  pref=paste(samps,collapse = "_vs_")
 system(paste0("mergePeaks -venn ",pref,".venn -matrix ",pref,".matrix -d given ",samps[1],"*homer.peaks ",samps[2],"*homer.peaks > ",pref,".mergedPeaks.txt"))
  
}

##### Run locally ####
#copy results to local folder
setwd("~/share/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/")
dir.create("merged_peaks_homer")
setwd("merged_peaks_homer/")
system("scp csuc:/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/merged_peaks/* .")

venn_files=list.files(pattern = "venn")
venn_files_diff_conditions=grep("vs",venn_files,value = T)
venn_files_diff_conditions=grep("ER",venn_files_diff_conditions,invert = T,value = T)

dir.create("~/share/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/merged_peaks_homer_venn_diagrams")

for (i in venn_files_diff_conditions) {
  venn=read.table(i,sep = "\t",header = T)
  compar=gsub(".venn","",i)
  name1=strsplit(compar,split = "_vs_")[[1]][1]
  name2=strsplit(compar,split = "_vs_")[[1]][2]
  
  pdf(paste0("~/share/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/merged_peaks_homer_venn_diagrams/",i,".pdf"))
  venn.plot <- draw.pairwise.venn(venn$Total[1]+venn$Total[3],
                                  venn$Total[2]+venn$Total[3], venn$Total[3], c(name1, name2));
  
  dev.off()
  
}

