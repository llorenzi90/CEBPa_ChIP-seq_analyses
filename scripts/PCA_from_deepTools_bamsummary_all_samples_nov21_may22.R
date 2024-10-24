## ---------------------------
##
##
## Purpose of script: pca of all ChIP-seq samples coveraged based
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-05-05
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
library(DESeq2)
library(ggplot2)
## ---------------------------
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/sample_correlation/all_samples_coverage_deepTools/")

da=read.table("whole_genome_readCounts.allbams.tab")

list_bf_whole_path=read.table("list_bam_files.txt")
list_bf=basename(list_bf_whole_path$V1)
length(unique(list_bf))
list_bf_whole_path$V1[duplicated(list_bf)]
#the duplicated ones are the old ones, the ones in which split failed so we 
#want to get rid of those, and it's easy to do so, just:
da_coords=da[,1:3]
da=da[,4:65]
da=da[,!duplicated(list_bf)]
list_bf=list_bf[!duplicated(list_bf)]
list_bf=gsub("sorted.markedDups.bam","nov21",list_bf)
list_bf=gsub("sorted.markedDups.proper_pairs.minq2.bam","apr22",list_bf)
colnames(da)=list_bf



libsizes <- colSums(da)
scaled_counts <- apply(da,1, function(x)x/libsizes*1000000)
scaled_counts <- scaled_counts[,colVars(scaled_counts)!=0]
pca <- prcomp(log2(scaled_counts+0.01),scale. = T)
dtp <- as.data.frame(pca$x)
dtp$samples <- rownames(dtp)
dtp$batch <- gsub("(.*)\\.(.*)","\\2",dtp$samples)
dtp$type <- "ChIP"
dtp$type[grepl("INP",dtp$samples)] <- "input"
dtp$vector <- ifelse(grepl("EV",dtp$samples,ignore.case = T),"EV",
                     ifelse(grepl("p42",dtp$samples,ignore.case = T),"p42","p30"))

dtp$factor <- ifelse(grepl("ER",dtp$samples,ignore.case = T),"ER",
                     ifelse(grepl("K4",dtp$samples,ignore.case = T),"K4","K27"))
dtp$LPS <- ifelse(grepl("LPS",dtp$samples,ignore.case = T),"LPS","UT")
dtp$replicate <- ifelse(grepl("R3",dtp$samples),"3",
                     ifelse(grepl("R2",dtp$samples),"2","1"))
percvar=round(pca$sdev^2/sum(pca$sdev^2)*100,2)
dir.create("PCA_plots")

dtp$factor_input=dtp$factor
dtp$factor_input[dtp$type=="input"]="input"

for (colvar in c("batch"    , "type", "vector"  ,  "factor",    "LPS" ,      "replicate","factor_input")) {
  
  pdf(paste0("PCA_plots/PCA_all_Samples_inclinput.colBy",colvar,".pdf"))
  g=ggplot(dtp,aes_string(x="PC1",y="PC2",col=colvar)) + 
    geom_point() +  geom_text(label=dtp$samples, nudge_x = 0.25, nudge_y = 0.25,
                              check_overlap = T) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
    ylab(paste0("PC2(",percvar[2],"% var)"))
  
  print(g)
  dev.off()
}

sampleDists <- dist(log2(scaled_counts+0.01))
#sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
library(RColorBrewer)
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#dir.create("heatmaps")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 8)
dir.create("heatmaps")
pdf("heatmaps/pheatmap_all_Samples_inclinput.log2+0.01.pdf")

print(p)
dev.off()

#separate data between each factor and do PCA

for (fac in unique(dtp$factor_input)) {
  tempdtp=dtp[dtp$factor_input==fac,]
  pdf(paste0("PCA_plots/PCA_all_Samples_",fac,".pdf"))
  g=ggplot(tempdtp,aes_string(x="PC1",y="PC2",col="vector")) + 
    geom_point(aes(shape=LPS)) +  geom_text(label=tempdtp$samples, nudge_x = 0.25, nudge_y = 0.25,
                              check_overlap = T) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
    ylab(paste0("PC2(",percvar[2],"% var)"))
  
  print(g)
  dev.off()
}
