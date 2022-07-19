setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/comparison_Jakobsen_et_al_2019_peaks/")
source("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/compute_plot_overlap.R")

peaks_paper=read.table("Jakobsen2019_CEBPa_peaks.mm10_liftOver.mm39_liftOver.bed")
peaks_our_p30=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/ER/peaks/nov21/p30_ER.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")
peaks_our_p42=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/ER/peaks/nov21/p42_ER.sorted.markedDups.bam.NoModel_q_0.05_peaks.narrowPeak")
peaks_our_merged=compute_plot_overlap(mylist = list(p30=peaks_our_p30,
                                                 p42=peaks_our_p42))[[1]]
peaks_paper_p30=peaks_paper[peaks_paper$V5%in%c("p30","common"),]
peaks_paper_p42=peaks_paper[peaks_paper$V5%in%c("p42","common"),]
peaks_paper_merged=compute_plot_overlap(mylist = list(p30=peaks_paper_p30,
                                                      p42=peaks_paper_p42))[[1]]
dir.create("venn_diagrams")
pdf("venn_diagrams/Overlap_Jakobsen2019_peaks.pdf",onefile = T)
overlap_all=compute_plot_overlap(mylist = list(our_peaks=peaks_our_merged,
                                               paper_peaks=peaks_paper_merged))
overlap_all[[2]]


####These include both common and p30/p42 peaks
overlap_p30=compute_plot_overlap(mylist = list(our_peaks=peaks_our_p30,
                                               paper_peaks=peaks_paper_p30),tit = "p30_our_peaks_vs_paper_peaks")
overlap_p30[[2]]


overlap_p42=compute_plot_overlap(mylist = list(our_peaks=peaks_our_p42,
                                               paper_peaks=peaks_paper_p42),tit = "p42_our_peaks_vs_paper_peaks")
overlap_p42[[2]]

#######comparison of p30 only and p42 only
overlap_p30_only=compute_plot_overlap(mylist = list(our_peaks=peaks_our_merged[peaks_our_merged$V6=="p30",],
                                                    paper_peaks=peaks_paper_merged[peaks_paper_merged$V6=="p30",]),tit = "p30_only_peaks_our_vs_paper")
overlap_p30_only[[2]]

overlap_p42_only=compute_plot_overlap(mylist = list(our_peaks=peaks_our_merged[peaks_our_merged$V6=="p42",],
                                                    paper_peaks=peaks_paper_merged[peaks_paper_merged$V6=="p42",]),tit = "p42_only_peaks_our_vs_paper")
overlap_p42_only[[2]]

######overlap all sets individually
overlap_all=compute_plot_overlap(mylist = list(p30_our=peaks_our_p30,
                                               p30_paper=peaks_paper_p30,
                                               p42_our=peaks_our_p42,
                                               p42_paper=peaks_paper_p42))
overlap_all[[2]]
dev.off()
