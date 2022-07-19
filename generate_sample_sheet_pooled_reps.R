setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/DiffBind/")
system("scp csuc:/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind/sample_sheet_csuc-narrow.csv .")

old_sample_sheet=read.csv("sample_sheet_csuc-narrow.csv",header = T)

#replace individual peak files by pooled reps files

system("ssh csuc wc -l /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/macs2_pooled_reps/*narrowPeak > ../nov21_pooled_peaks_numbers.txt")
system("ssh csuc wc -l /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/macs2_pooled_reps/may22/*narrowPeak > ../may22_pooled_peaks_numbers.txt")
system("ssh csuc wc -l /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/macs2_pooled_reps/may22/*narrowPeak > ../may22_pooled_peaks_numbers.txt")


nov21_pooled=read.table("../nov21_pooled_peaks_numbers.txt")
nov21_pooled=nov21_pooled[-nrow(nov21_pooled),]
nov21_pooled$sample_id=gsub(".pooled.NoModel_q_0.05_peaks.narrowPeak","",basename(nov21_pooled$V2))
old_sample_sheet$id=paste0(old_sample_sheet$Condition,"_",old_sample_sheet$Treatment)
table(nov21_pooled$sample_id%in%old_sample_sheet$id)
nov21_pooled$sample_id[!nov21_pooled$sample_id%in%old_sample_sheet$id]

old_sample_sheet$Peaks=nov21_pooled$V2[match(old_sample_sheet$id,
                                             nov21_pooled$sample_id)]
old_sample_sheet$Replicate=gsub("R","",old_sample_sheet$Replicate)

write.csv(old_sample_sheet[,-ncol(old_sample_sheet)],"sample_sheet_csuc.pooledReps.csv",row.names = F)
