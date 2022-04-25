# y la otra cosa es: de los peaks que salen up y down en la comparación p30LPSvsp42LPS, se le podrían restar los picos que ya hayan salido up y down en p30UTvsp42UT?
#   el motivo es que muchos son los mismos, y nos gustaría ver específicamente los que salen sólo en la comparación LPS
comp1 <- "p30LPSvsp42LPS"
comp2 <- "p30UTvsp42UT"

comps=c(comp1,comp2)

deseq2dir <- "~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/DESeq2_results/"
list.files(deseq2dir)
deseq2res <- lapply(comps, function(x){
  read.csv(paste0(deseq2dir,x,"_DESeq2_res.noNAs.csv"))
})

names(deseq2res) <-comps
padjco=0.05
diffs <- lapply(deseq2res, function(x){
  ifelse(x$padj<=padjco,ifelse(x$log2FoldChange>0,"UP","DOWN"),"NO")
})
table(diffs$p30LPSvsp42LPS)

downs <- lapply(names(diffs), function(x){
  deseq2res[[x]]$X[diffs[[x]]=="DOWN"]
})

ups <- lapply(names(diffs), function(x){
  deseq2res[[x]]$X[diffs[[x]]=="UP"]
})

library(ggvenn)
ggvenn(data = downs)
library(RVenn)
library(VennDiagram)
names(ups) <- comps
names(downs) <- comps
ggvenn::ggvenn(ups) + ggtitle("UP peaks")
ggvenn::ggvenn(downs) + ggtitle("DOWN peaks")

up_LPS_only <- ups$p30LPSvsp42LPS[!ups$p30LPSvsp42LPS%in%ups$p30UTvsp42UT]
down_LPS_only <- downs$p30LPSvsp42LPS[!downs$p30LPSvsp42LPS%in%downs$p30UTvsp42UT]

peaks_file <- "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/macs2_pooled_reps/merged_peaks/merged_peaks_EV.p30.p42_LPS.UT.sorted.bed.withcolnames.annotated.csv"
peaks <- read.csv(peaks_file)
cols_to_keep <-c(1:3,19,8,5) 
head(peaks[,cols_to_keep])
peaks$strand <- "."

backgroundPeaks <- deseq2res$p30LPSvsp42LPS$X[!deseq2res$p30LPSvsp42LPS$X%in%c(up_LPS_only,down_LPS_only)]
list_peak.sets <- list(upPeaks_LPSonly=up_LPS_only,
                       downPeaks_LPSonly=down_LPS_only,
                       backgroundPeaks=backgroundPeaks)

datadir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/DESeq2_results/"

setwd(datadir)
comp="p30LPSvsp42LPS_LPSonly_Peaks"
dir.create(comp)
setwd(comp)


pdf("Up_peaks_overlap_LPS_UT.pdf")
ggvenn::ggvenn(ups) + ggtitle("UP peaks")
dev.off()

pdf("Down_peaks_overlap_LPS_UT.pdf")
ggvenn::ggvenn(downs) + ggtitle("DOWN peaks")
dev.off()

outdir="bed_files_homer_input/"
dir.create(outdir)
for (peakset in names(list_peak.sets)) {
  peakids <- list_peak.sets[[peakset]]
  
  write.table(peaks[peaks$peakID%in%peakids,cols_to_keep],
              paste0(outdir,comp,"_",peakset,".bed"),
              col.names = F,row.names = F,quote = F, sep = "\t")
}
