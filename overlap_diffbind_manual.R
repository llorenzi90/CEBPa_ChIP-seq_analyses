## ---------------------------
##
##
## Purpose of script:comparison DiffBind peaks
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-03-31
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


require(DiffBind)
## ---------------------------
setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/DiffBind/DESeq2/")
Dba3first=readRDS("Dba_first_3comps.RDS")

diffbind_peaks=dba.peakset(Dba3first,bRetrieve = T,
                           DataType = DBA_DATA_FRAME)

diffbind_peaks$ID=paste0(diffbind_peaks$CHR,":",diffbind_peaks$START,"-",diffbind_peaks$END)
manual_peaks=read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/macs2_pooled_reps/merged_peaks/merged_peaks_EV.p30.p42_LPS.UT.sorted.bed.withcolnames.annotated.csv")
manual_peaks$peakID
table(diffbind_peaks$ID%in%manual_peaks$peakID)
diffbind_peaks$width=diffbind_peaks$END - diffbind_peaks$START +1
summary(diffbind_peaks$width)
summary(manual_peaks$width)
plot(density(diffbind_peaks$width))
points(density(manual_peaks$width),type="s",col="red")

plot(density(manual_peaks$width))
points(density(diffbind_peaks$width),type="s",col="red")

plot(density(log10(diffbind_peaks$width)))
points(density(log10(manual_peaks$width)),type = "s",col="darkred")
library(bedtoolsr)
diffbind_peaks$START=diffbind_peaks$START-1
manual_peaks$start=manual_peaks$start-1
inters=bt.intersect(diffbind_peaks[,c(1:3,16:17)],manual_peaks,wao = T)
colnames(inters)=c(paste0("DB_",colnames(diffbind_peaks)[c(1:3,16:17)]),
                   paste0("M_",colnames(manual_peaks)),
                   "overlap")

peaks_with_overlap=inters[inters$overlap!=0,]
library(gplots)

plot(log2(as.integer(peaks_with_overlap$M_width)),log2(peaks_with_overlap$DB_width))
#
tDB=table(peaks_with_overlap$DB_ID)
table(tDB)
table(tDB>1)
table(tDB>1)/nrow(diffbind_peaks)*100
table(tDB>1)/length(tDB)*100
tM=table(peaks_with_overlap$M_peakID)
table(tM)
table(tM>1)
table(tM>1)/nrow(manual_peaks)*100
length(unique(peaks_with_overlap$M_peakID))/nrow(manual_peaks)*100
table(tM>1)/length(tM)*100

tM[tM==9]
#compare length
#compare matching peaks X overlap
#check how many diffbind peaks contain many manual peaks and viceversa

library(ggvenn)
Dba3first$contrasts[[1]]$contrast
Dba3first$DESeq2
resDB=as.data.frame(dba.report(Dba3first,contrast = 1,th = 1))
resDB$ID=paste0(resDB$seqnames,":",resDB$start,"-",resDB$end)
table(resDB$ID%in%diffbind_peaks$ID)
DBpeaksDB=resDB$ID[resDB$FDR<=0.05]
table(DBpeaksDB%in%peaks_with_overlap$DB_ID)

resM=read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/DESeq2_results/annotated_results/p30UTvsp42UT_DESeq2_res.noNAs.annotated.csv")
table(resM$X%in%manual_peaks$peakID)
resM$X[!resM$X%in%manual_peaks$peakID]
DBpeaksM=resM$X[resM$padj<=0.05]
table(DBpeaksM%in%manual_peaks$peakID)
resM=resM[resM$X%in%manual_peaks$peakID,]
DBpeaksM=resM$X[resM$padj<=0.05]

table(DBpeaksM%in%peaks_with_overlap$M_peakID)

length(unique(peaks_with_overlap$DB_ID[peaks_with_overlap$M_peakID%in%DBpeaksM]))
length(unique(peaks_with_overlap$M_peakID[peaks_with_overlap$DB_ID%in%DBpeaksDB]))

diff_peaks_with_overlapM=peaks_with_overlap[peaks_with_overlap$M_peakID%in%DBpeaksM,]
table(table(diff_peaks_with_overlapM$DB_ID))

diff_peaks_with_overlapM=peaks_with_overlap[peaks_with_overlap$M_peakID%in%DBpeaksM,]
table(table(diff_peaks_with_overlapM$DB_ID))
table(table(diff_peaks_with_overlapM$M_peakID))
table(table(diff_peaks_with_overlapM$M_peakID)>1)


diff_peaks_with_overlapDB=peaks_with_overlap[peaks_with_overlap$DB_ID%in%DBpeaksDB,]
table(table(diff_peaks_with_overlapDB$DB_ID))
table(table(diff_peaks_with_overlapDB$M_peakID))
table(table(diff_peaks_with_overlapM$M_peakID)>1)

table(unique(diff_peaks_with_overlapM$M_peakID)%in%diff_peaks_with_overlapDB$M_peakID)
table(unique(diff_peaks_with_overlapM$M_peakID)%in%diff_peaks_with_overlapDB$M_peakID)

table(unique(diff_peaks_with_overlapDB$DB_ID)%in%diff_peaks_with_overlapM$DB_ID)

up_peaksDB=resDB$ID[resDB$ID%in%DBpeaksDB&resDB$Fold>0]
down_peaksDB=resDB$ID[resDB$ID%in%DBpeaksDB&resDB$Fold<0]

up_peaksM=resM$X[resM$X%in%DBpeaksM&resM$log2FoldChange>0]
down_peaksM=resM$X[resM$X%in%DBpeaksM&resM$log2FoldChange<0]

up_down_lists=list(upPeaks=list(DB=up_peaksDB,M=up_peaksM),
                   downPeaks=list(DB=down_peaksDB,M=down_peaksM))

for (na in names(up_down_lists)) {
  print("peakset: ")
  print(na)
  ps=up_down_lists[[na]]
  print(paste("DiffBind", na,"with overlap manual peaks"))
  print(table(ps$DB%in%peaks_with_overlap$DB_ID))
  print(table(ps$DB%in%peaks_with_overlap$DB_ID)/length(ps$DB)*100)
  print(paste("Manual", na,"with overlap DiffBind peaks"))
  print(table(ps$M%in%peaks_with_overlap$M_peakID))
  print(table(ps$M%in%peaks_with_overlap$M_peakID)/length(ps$M)*100)
  DB_peaksWO=peaks_with_overlap[peaks_with_overlap$DB_ID%in%ps$DB,]
  M_peaksWO=peaks_with_overlap[peaks_with_overlap$M_peakID%in%ps$M,]
  print(paste("DiffBind", na,"with overlap manual diff peaks"))
  print(table(ps$DB%in%diff_peaks_with_overlapM$DB_ID))
  print(table(ps$DB%in%diff_peaks_with_overlapM$DB_ID)/length(ps$DB)*100)
  print(paste("DiffBind", na,"with overlap manual",na))
  print(table(ps$DB%in%M_peaksWO$DB_ID))
  print(table(ps$DB%in%M_peaksWO$DB_ID)/length(ps$DB)*100)
  print(paste("Manual", na,"with overlap DiffBind diff peaks"))
  print(table(ps$M%in%diff_peaks_with_overlapDB$M_peakID))
  print(table(ps$M%in%diff_peaks_with_overlapDB$M_peakID)/length(ps$M)*100)
  print(paste("Manual", na,"with overlap DiffBind",na))
  print(table(ps$M%in%DB_peaksWO$M_peakID))
  print(table(ps$M%in%DB_peaksWO$M_peakID)/length(ps$M)*100)
  print(paste("DiffBind",na,"spanning more than 1 manual peak"))
  print(table(table(DB_peaksWO$DB_ID)>1))
  print(table(table(DB_peaksWO$DB_ID)>1)/length(ps$DB)*100)
  print(paste("Manual",na,"spanning more than 1 DiffBind peak"))
  print(table(table(M_peaksWO$M_peakID)>1))
  print(table(table(M_peaksWO$M_peakID)>1)/length(ps$M)*100)
 
}

gl=readRDS("Dba_remainingcomps.greylist.RDS")
diffbind_peaksGL=dba.peakset(gl,bRetrieve = T,
                           DataType = DBA_DATA_FRAME)

diffbind_peaksGL$ID=paste0(diffbind_peaksGL$CHR,":",diffbind_peaksGL$START,"-",diffbind_peaksGL$END)
blpeaks=diffbind_peaks$ID[!diffbind_peaks$ID%in%diffbind_peaksGL$ID]
write.table(blpeaks,"/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/DiffBind/greylisted_peaks.txt",row.names = F,col.names = F,quote = F)
write.table(peaks_with_overlap$M_peakID[peaks_with_overlap$DB_ID%in%blpeaks],"/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/DiffBind/greylisted_peaks_overlap_manualapproach.peakID.txt",row.names = F,col.names = F,quote = F)

blpeaks_M=peaks_with_overlap$M_peakID[peaks_with_overlap$DB_ID%in%blpeaks]



table(blpeaks%in%peaks_with_overlap$DB_ID)
table(blpeaks%in%diff_peaks_with_overlapM$DB_ID)


View(diffbind_peaks[diffbind_peaks$ID%in%blpeaks,c(4:15)])
diffbind_peaks_counts=diffbind_peaks[,4:15]
libs=colSums(diffbind_peaks_counts)
cpmdiffbind=apply(diffbind_peaks_counts, 1,function(x)x/libs*1000000)
cpmdiffbind=t(cpmdiffbind)
head(cpmdiffbind)
View(cpmdiffbind[diffbind_peaks$ID%in%blpeaks,])

manual_counts=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/merged_peaks_EV.p30.p42_LPS.UT.counts")
libs_manual=colSums(manual_counts)
libs/1000000
libs_manual/1000000
cpm_manual=t(apply(manual_counts, 1, function(x)x/libs_manual*1000000))
View(cpm_manual[rownames(cpm_manual)%in%blpeaks_M,])
View(cpm_manual)
#Compare greylisted with both diffbind and manual
#compare number of differential peaks
#how many differential peaks are unique for each approach?
#how many differential peaks are marked as greylist?
#check if diffbind provides a way of checking consistency between replicates (dba.plotVenn)
#read diffbind manual!
olap.rate <- dba.overlap(Dba3first,mode=DBA_OLAP_RATE)
plot(olap.rate,type='b',ylab='# peaks',
     xlab='Overlap at least this many peaksets')

#How well do the replicates agree on their peak calls?
Dba3first$masks$p30&Dba3first$masks$LPS
dba.overlap(Dba3first,Dba3first$masks$p30&Dba3first$masks$LPS,
            mode=DBA_OLAP_RATE)
dev.new()
dba.plotVenn(Dba3first,Dba3first$masks$p30&Dba3first$masks$LPS)
dba.plotVenn(Dba3first,Dba3first$masks$p30&Dba3first$masks$UT)
dba.plotVenn(Dba3first,Dba3first$masks$p42&Dba3first$masks$LPS)
dba.plotVenn(Dba3first,Dba3first$masks$p42&Dba3first$masks$UT)
dba.plotVenn(Dba3first,c(6,7))
colSums(diffbind_peaks[,4:15])/1000000

Dba3first_consensus <- dba.peakset(Dba3first,
                                     consensus=DBA_CONDITION,
                                     minOverlap=0.66)
#what are each type of duplicates? optical, nonoptical?

#Things to check:
#why diffbind counts are so different?
#I know now: I was using bUseSummarizeOverlaps=F, this has to be set to TRUE if we want it
# to consider paired end reads. see https://support.bioconductor.org/p/58135/ I need to check all this parameters in more detail
#normalization
#consensus peaks

#If I have time left: try to do clustering based on peaks called only

