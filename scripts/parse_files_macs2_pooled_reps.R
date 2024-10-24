samps=read.table("~/share/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/scripts/sample_control_macs2_Maria.txt")
library(tidyr)
samps1=separate(samps,col = 1,into = c(NA,NA,NA,NA,"sample","bam"),sep = "/",remove = F)
samps2=separate(samps,col = 2,into = c(NA,NA,NA,NA,"sample","bam"),sep = "/",remove = F)
colnames(samps1)[1]="path"
colnames(samps2)[2]="path"
samps=rbind(samps1[,1:3],samps2[,2:4])
#update bam files
newbams=read.table("~/share/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/scripts/new_bam_files_Maria.txt")
newbams <- separate(newbams,col = 1,into = c(NA,NA,NA,NA,NA,"sample","bam"),sep = "/",remove = F)

table(newbams$sample%in%samps$sample)
samps$path[match(newbams$sample,samps$sample)] <- newbams$V1

#remove rep from sample name

samps$sample <- gsub("P42","p42",samps$sample)
samps$repl=gsub("^_","",gsub("(.*)(_R[1-9]*)","\\2",samps$sample))
samps$repl[grep("INP",samps$sample)] <- paste0("INP_",samps$repl[grep("INP",samps$sample)])
samps$repl

samps$csamp <- gsub("_R[1-9]*","",gsub("INP_","",samps$sample))
samps$repl

#split sampls in input and not input
inps=samps[grep("INP",samps$sample),]
samps=samps[grep("INP",samps$sample,invert = T),]
table(samps$csamp%in%inps$csamp)

reps=c()
inpreps=c()
for (sa in unique(samps$csamp)) {
  reps=c(reps,paste0(samps$path[samps$csamp==sa],collapse = ","))
  inpreps=c(inpreps,paste0(inps$path[inps$csamp==sa],collapse = ","))
}

df_out=data.frame(sample=unique(samps$csamp),treat=reps,control=inpreps)
write.table(df_out,"~/share/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/scripts/macs2_pooled_sample_sheet.txt",quote = F,col.names = F,row.names = F,sep = " ")
