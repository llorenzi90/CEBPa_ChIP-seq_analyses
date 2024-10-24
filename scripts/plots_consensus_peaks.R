manual_peaks=read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/macs2_pooled_reps/merged_peaks/merged_peaks_EV.p30.p42_LPS.UT.sorted.bed.withcolnames.annotated.csv")
par(las=1)
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/")
dir.create("consensus_peaks_plots")

pdf("consensus_peaks_plots/Nsamples.vs.Npeaks.pdf")
par(las=1)
barplot(table(manual_peaks$nsamples)/1000,xlab = "Peaks found in as many samples",
        cex.names = 2,cex.axis = 2,cex.lab=2)
dev.off()

pdf("consensus_peaks_plots/Nsamples.vs.q-val.pdf")
par(las=1)
boxplot(split(x = 10^(-manual_peaks$meanscore),f = manual_peaks$nsamples),
        xlab="Peaks found in as many samples",cex.names = 2,cex.axis = 2,cex.lab=2)
dev.off()

manual_peaks$simplified_annot <- manual_peaks$annotation
manual_peaks$simplified_annot[grep("Exon",manual_peaks$annotation)] <- "exon"
manual_peaks$simplified_annot[grepl("Exon",manual_peaks$annotation)&
                                     grepl("exon 1 of",manual_peaks$annotation)] <- "exon (1st)"

manual_peaks$simplified_annot[grep("Intron",manual_peaks$annotation)] <- "intron"

manual_peaks$simplified_annot[grepl("Intron",manual_peaks$annotation)&
                                     grepl("intron 1 of",manual_peaks$annotation)] <- "intron (1st)"
table(manual_peaks$simplified_annot)
topie <- table(manual_peaks$simplified_annot)
source("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/cbPalette.R")
cols=colorRampPalette( cbPalette)(length(topie))
names(cols)=names(topie)
library(ggplot2)
library(ggrepel)
library(tidyverse)

df=data.frame(sort(topie))

df2 <- df %>% 
  mutate(csum = rev(cumsum(rev(Freq))), 
         pos = Freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Freq/2, pos))

g=ggplot(df, aes(x = "" , y = Freq, fill = fct_inorder(Var1))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = cols) +
  geom_label_repel(data = df2,
                   aes(y = pos, label = Var1),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  theme_void() +theme(legend.position = "none") 
pdf("consensus_peaks_plots/Consensus_peaks_anntotation.pie.pdf")       
print(g)
dev.off()

pdf("consensus_peaks_plots/peak_length_distribution.pdf")
ggplot(manual_peaks,aes(x=width))+geom_density()+
  scale_x_log10()+theme_classic()+
  xlab("peak length (bp)")
dev.off()

