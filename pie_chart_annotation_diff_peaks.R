## Notes: First version of this script was WRONG!!! I used the human annotation instead of mouse annotation!!
##   This has been corrected now (20 may 22)

p30UTvsp42UT_DESeq2_res.noNAs.annotated <- read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/DESeq2_results/annotated_results/p30UTvsp42UT_DESeq2_res.noNAs.annotated.csv")
p30UTvsp42UT_DESeq2_res.noNAs.annotated <- read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/K27/macs2_merged_peaks_counts/DESeq2_results/annotated_results/p30UTvsp42UT_DESeq2_res.noNAs.annotated.csv")

p30UTvsp42UT_DESeq2_res.noNAs.annotated$simplified_annot <- p30UTvsp42UT_DESeq2_res.noNAs.annotated$annotation
p30UTvsp42UT_DESeq2_res.noNAs.annotated$simplified_annot[grep("Exon",p30UTvsp42UT_DESeq2_res.noNAs.annotated$annotation)] <- "exon"
p30UTvsp42UT_DESeq2_res.noNAs.annotated$simplified_annot[grepl("Exon",p30UTvsp42UT_DESeq2_res.noNAs.annotated$annotation)&
                                grepl("exon 1 of",p30UTvsp42UT_DESeq2_res.noNAs.annotated$annotation)] <- "exon (1st)"

p30UTvsp42UT_DESeq2_res.noNAs.annotated$simplified_annot[grep("Intron",p30UTvsp42UT_DESeq2_res.noNAs.annotated$annotation)] <- "intron"

p30UTvsp42UT_DESeq2_res.noNAs.annotated$simplified_annot[grepl("Intron",p30UTvsp42UT_DESeq2_res.noNAs.annotated$annotation)&
                                grepl("intron 1 of",p30UTvsp42UT_DESeq2_res.noNAs.annotated$annotation)] <- "intron (1st)"
table(p30UTvsp42UT_DESeq2_res.noNAs.annotated$simplified_annot)

p30UTvsp42UT_DESeq2_res.noNAs.annotated$diff <- ifelse(p30UTvsp42UT_DESeq2_res.noNAs.annotated$padj<=0.05,ifelse(p30UTvsp42UT_DESeq2_res.noNAs.annotated$log2FoldChange>0,"UP","DOWN"),"NO")

table(p30UTvsp42UT_DESeq2_res.noNAs.annotated$diff,p30UTvsp42UT_DESeq2_res.noNAs.annotated$simplified_annot!="Promoter (<=1kb)")
table(p30UTvsp42UT_DESeq2_res.noNAs.annotated$simplified_annot,p30UTvsp42UT_DESeq2_res.noNAs.annotated$diff)
source("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/cbPalette.R")
library(ggplot2)
library(ggrepel)
library(tidyverse)
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/K27/macs2_merged_peaks_counts/DESeq2_results/")
dir.create("pie_charts_annotation_differential_peaks")
setwd("pie_charts_annotation_differential_peaks/")
dir.create("p30UTvsp42UT")
setwd("p30UTvsp42UT/")
for (set in c("UP","NO","DOWN")) {
  topie <- table(p30UTvsp42UT_DESeq2_res.noNAs.annotated$simplified_annot[p30UTvsp42UT_DESeq2_res.noNAs.annotated$diff==set])
  print(topie/sum(topie)*100)
  cols=colorRampPalette( cbPalette)(length(topie))
  names(cols)=names(topie)
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
  
  pdf(paste0(set,"diffbind_anntotation.pie.pdf") )      
  print(g)
  dev.off()
  

}




