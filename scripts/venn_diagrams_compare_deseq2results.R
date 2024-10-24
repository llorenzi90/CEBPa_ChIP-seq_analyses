## ---------------------------
##
##
## Purpose of script:compare peaks up in p30vsEV and in p42vsEV
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-04-03
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
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/DESeq2_results/annotated_results/")

p30vsEV=read.csv("p30UTvsEVUT_DESeq2_res.noNAs.annotated.csv")
p42vsEV=read.csv("p42UTvsEVUT_DESeq2_res.noNAs.annotated.csv")
padjco=0.05
reslist=list(p30vsEV=p30vsEV,p42vsEV=p42vsEV)
up_peaks=lapply(reslist, function(x)x$X[x$padj<=padjco&x$log2FoldChange>0])
down_peaks=lapply(reslist, function(x)x$X[x$padj<=padjco&x$log2FoldChange<0])

library(gplots)
v.up <- venn(up_peaks,show.plot = F) #this function also makes a venn diagram (but an ugly one)
#but I indicate not to show the plots and store the results in v.up object
v.down <- venn(down_peaks,show.plot = F)

source("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/cbPalette.R")
library(ggvenn)
setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/DESeq2_results/")
dir.create("venn_diagrams")
setwd("venn_diagrams/")
dir.create("comparison_p30vsEV_vs_p42vsEV")
setwd("comparison_p30vsEV_vs_p42vsEV/")
for (si in c(6,7,8,9,10)) {
  g <- ggvenn(
    up_peaks, 
    fill_color = cbPalette[1:2], #change this to set other colors 
    #fill_color = c( "#EFC000FF", "#868686FF"),
    stroke_size = 0.5, set_name_size = si,
    show_percentage = F,
    text_size = si
  ) + ggtitle("Up peaks") + theme(plot.title =  element_text(size=si*2.5))
  pdf(paste0("upPeaks_text_size_",si,".pdf"))
  print(g)
  dev.off()
  
  g <- ggvenn(
    down_peaks, 
    fill_color = cbPalette[1:2], #change this to set other colors 
    #fill_color = c( "#EFC000FF", "#868686FF"),
    stroke_size = 0.5, set_name_size = si,
    show_percentage = F,
    text_size = si
  ) + ggtitle("Down peaks") + theme(plot.title =  element_text(size=si*2.5))
  pdf(paste0("downPeaks_text_size_",si,".pdf"))
  print(g)
  dev.off()
}


setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/DESeq2_results/")

ress=list.files("annotated_results/")
combinations=combn(ress,m=2)
combinations[,1]
#Do all bi-comparisons:
for (c in 1:ncol(combinations)) {
  comp=combinations[,c]
  reslist=lapply(comp, function(x)read.csv(paste0("annotated_results/",x)))
  names(reslist) <- sapply(comp,function(x)gsub("_DESeq2_res.noNAs.annotated.csv","",x))
  comp=paste(names(reslist),collapse = "_vs_")
  outdir=paste0("venn_diagrams/",comp)
  dir.create(outdir) 
  up_peaks=lapply(reslist, function(x)x$X[x$padj<=padjco&x$log2FoldChange>0])
  down_peaks=lapply(reslist, function(x)x$X[x$padj<=padjco&x$log2FoldChange<0])
  
  for (si in c(6,7,8,9,10)) {
    g <- ggvenn(
      up_peaks, 
      fill_color = cbPalette[1:2], #change this to set other colors 
      #fill_color = c( "#EFC000FF", "#868686FF"),
      stroke_size = 0.5, set_name_size = si,
      show_percentage = F,
      text_size = si
    ) + ggtitle("Up peaks") + theme(plot.title =  element_text(size=si*2.5))
    pdf(paste0(outdir,"/upPeaks_text_size_",si,".pdf"))
    print(g)
    dev.off()
    
    g <- ggvenn(
      down_peaks, 
      fill_color = cbPalette[1:2], #change this to set other colors 
      #fill_color = c( "#EFC000FF", "#868686FF"),
      stroke_size = 0.5, set_name_size = si,
      show_percentage = F,
      text_size = si
    ) + ggtitle("Down peaks") + theme(plot.title =  element_text(size=si*2.5))
    pdf(paste0(outdir,"/downPeaks_text_size_",si,".pdf"))
    print(g)
    dev.off()
  }
  
  #add extra comparison:
  reslist[[2]]$log2FoldChange <- - reslist[[2]]$log2FoldChange
  names(reslist)[2] <- paste(strsplit(names(reslist)[2],split = "vs")[[1]][2:1],collapse = "vs")
  comp=paste(names(reslist),collapse = "_vs_")
  outdir=paste0("venn_diagrams/",comp)
  dir.create(outdir) 
  up_peaks=lapply(reslist, function(x)x$X[x$padj<=padjco&x$log2FoldChange>0])
  down_peaks=lapply(reslist, function(x)x$X[x$padj<=padjco&x$log2FoldChange<0])
  for (si in c(6,7,8,9,10)) {
    g <- ggvenn(
      up_peaks, 
      fill_color = cbPalette[1:2], #change this to set other colors 
      #fill_color = c( "#EFC000FF", "#868686FF"),
      stroke_size = 0.5, set_name_size = si,
      show_percentage = F,
      text_size = si
    ) + ggtitle("Up peaks") + theme(plot.title =  element_text(size=si*2.5))
    pdf(paste0(outdir,"/upPeaks_text_size_",si,".pdf"))
    print(g)
    dev.off()
    
    g <- ggvenn(
      down_peaks, 
      fill_color = cbPalette[1:2], #change this to set other colors 
      #fill_color = c( "#EFC000FF", "#868686FF"),
      stroke_size = 0.5, set_name_size = si,
      show_percentage = F,
      text_size = si
    ) + ggtitle("Down peaks") + theme(plot.title =  element_text(size=si*2.5))
    pdf(paste0(outdir,"/downPeaks_text_size_",si,".pdf"))
    print(g)
    dev.off()
  }
}


#special comparison p30UTvsp42UT vs p42