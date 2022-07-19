setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/")
dir.create("macs2_may22")
setwd("macs2_may22/")
system("ssh csuc ls /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind/new_bam_files_may22/*bam > list_bam_files_CEBPa_may22.txt")
bam_files=read.table("list_bam_files_CEBPa_may22.txt")
bam_files$sampleID=gsub(".sorted.markedDups.proper_pairs.minq2.bam","",basename(bam_files$V1))
bam_files$norepnoinpID=gsub("INP_|R[1-4]_","",bam_files$sampleID)
length(unique(bam_files$norepnoinpID))
bam_files$norepnoinpID[grepl("INP",bam_files$V1)]%in%bam_files$norepnoinpID[!grepl("INP",bam_files$V1)]

bam_files_samples=bam_files[!grepl("INP",bam_files$V1),]
bam_files_input=bam_files[grepl("INP",bam_files$V1),]
bam_files_samples$input=bam_files_input$V1[match(bam_files_samples$norepnoinpID,
                                              bam_files_input$norepnoinpID)]

write.table(bam_files_samples[,c("V1","input")],"sample_control_CEBPa_may22.txt",row.names = F,quote = F,col.names = F,sep = " ")

#submit in cluster:
#cat sample_control_CEBPa_may22.txt |while read s c; do echo "sample $s"; echo -e "control $c";bn=$(basename $s);bn=${bn/.sorted.markedDups.proper_pairs.minq2.bam/};echo $bn;sbatch -J $bn.macs2 macs2_NoModel_q_0.05.sh $s $c ;done

#make file for pooled reps:
#format sampleID path1,path1 inputpath1,inputpath2

#wdir=/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/macs2_pooled_reps/may22

#cat pooled_sample_sheet_CEBPa_may22.txt |while read sample treat input;do echo $sample;echo $treat;echo $input;sbatch -J $sample.macs2_pooled run_macs2_pooled.sh $wdir $sample $treat $input;done

commasamples=aggregate(bam_files_samples$V1,by=list(sampleID=bam_files_samples$norepnoinpID),function(x)paste0(x,collapse = ","))
commainputs=aggregate(bam_files_input$V1,by=list(sampleID=bam_files_input$norepnoinpID),function(x)paste0(x,collapse = ","))
commasamples$input=commainputs$x[match(commasamples$sampleID,commainputs$sampleID)]
write.table(commasamples,"pooled_sample_sheet_CEBPa_may22.txt",quote = F,col.names = F,row.names = F,sep = " ")
