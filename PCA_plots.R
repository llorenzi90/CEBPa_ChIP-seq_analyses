## ---------------------------
##
##
## Purpose of script:DESeq2 CEBEPa
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-02-28 - Updated 03/04/2022
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
library(DESeq2)
library(ggplot2)
countdata=read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/merged_peaks_EV.p30.p42_LPS.UT.counts")

colnames(countdata) <- gsub(".sorted.markedDups.bam","",colnames(countdata))
colnames(countdata)

#quick PCA
libsizes <- colSums(countdata)
scaled_counts <- apply(countdata,1, function(x)x/libsizes*1000000)
pca <- prcomp(scaled_counts,scale. = T)
dtp <- as.data.frame(pca$x)
dtp$samples <- rownames(dtp)
dtp$vector=gsub("(.*)_(.*)_(.*)","\\1",dtp$samples)
dtp$vector[dtp$vector=="P42"] <- "p42"
dtp$LPS=gsub("(.*)_(.*)_(.*)","\\2",dtp$samples)
dtp$rep=gsub("(.*)_(.*)_(.*)","\\3",dtp$samples)

#percent variance:
percvar=round(pca$sdev^2/sum(pca$sdev^2)*100,2)

setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/")
dir.create("PCA_plots")

pdf("PCA_plots/PCA_CPM_merged_peaks_EV_p30_p42_LPS_UT.pdf",onefile = T)
ggplot(dtp, aes(x=PC1,y=PC2,col=vector)) + 
  geom_point() +theme_classic() + xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
ggplot(dtp, aes(x=PC1,y=PC2,col=LPS)) + 
  geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
ggplot(dtp, aes(x=PC1,y=PC2,col=rep)) + 
  geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
dev.off()
pdf("PCA_plots/PCA_CPM_merged_peaks_EV_p30_p42_LPS_UT.colLPS.pdf",onefile = T)
ggplot(dtp, aes(x=PC1,y=PC2,col=LPS)) + 
  geom_point(aes(shape=vector),size=4) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))+theme(text = element_text(size = 20))
dev.off()
##Add pca using 500 most variable peaks (this is what DESeq2 package does by default and
#what I have been doing for other analyses) as a comparison
rv <- rowVars(t(scaled_counts))
ntop=500 #select only the top 500 genes with highest variance
#This is optional, is the default of the DESeq2 package
select <- order(rv, decreasing = TRUE)[1:ntop]
top500variable <- scaled_counts[,select]
pca <- prcomp(top500variable,scale. = T)
dtp <- as.data.frame(pca$x)
dtp$samples <- rownames(dtp)
dtp$vector=gsub("(.*)_(.*)_(.*)","\\1",dtp$samples)
dtp$vector[dtp$vector=="P42"] <- "p42"
dtp$LPS=gsub("(.*)_(.*)_(.*)","\\2",dtp$samples)
dtp$rep=gsub("(.*)_(.*)_(.*)","\\3",dtp$samples)

#percent variance:
percvar=round(pca$sdev^2/sum(pca$sdev^2)*100,2)

pdf("PCA_plots/PCA_CPM_merged_peaks_EV_p30_p42_LPS_UT.top500moreVarPeaks.pdf",onefile = T)
ggplot(dtp, aes(x=PC1,y=PC2,col=vector)) + 
  geom_point() +theme_classic() + xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
ggplot(dtp, aes(x=PC1,y=PC2,col=LPS)) + 
  geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
ggplot(dtp, aes(x=PC1,y=PC2,col=rep)) + 
  geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
dev.off()

#However 500 peaks may be too few for the amount of peaks we have, try with % of peaks
fraction_of_peaks=c(0.01,0.05,0.10,0.20,0.30,0.40,0.50,0.60,0.70)

getwd()
dir.create("PCA_plots/varying_fraction_top_peaks")
rv <- rowVars(t(scaled_counts))
for(fr in fraction_of_peaks){
  ntop=round(ncol(scaled_counts)*fr) #select only the top 500 genes with highest variance
  #This is optional, is the default of the DESeq2 package
  select <- order(rv, decreasing = TRUE)[1:ntop]
  topXvariable <- scaled_counts[,select]
  pca <- prcomp(topXvariable,scale. = T)
  dtp <- as.data.frame(pca$x)
  dtp$samples <- rownames(dtp)
  dtp$vector=gsub("(.*)_(.*)_(.*)","\\1",dtp$samples)
  dtp$vector[dtp$vector=="P42"] <- "p42"
  dtp$LPS=gsub("(.*)_(.*)_(.*)","\\2",dtp$samples)
  dtp$rep=gsub("(.*)_(.*)_(.*)","\\3",dtp$samples)
  #percent variance:
  percvar=round(pca$sdev^2/sum(pca$sdev^2)*100,2)
  
  #generate ggplots
  g1=ggplot(dtp, aes(x=PC1,y=PC2,col=vector)) + 
    geom_point() +theme_classic() + xlab(paste0("PC1(",percvar[1],"% var)"))+
    ylab(paste0("PC2(",percvar[2],"% var)"))
  g2=ggplot(dtp, aes(x=PC1,y=PC2,col=LPS)) + 
    geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
    ylab(paste0("PC2(",percvar[2],"% var)"))
  g3=ggplot(dtp, aes(x=PC1,y=PC2,col=rep)) + 
    geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
    ylab(paste0("PC2(",percvar[2],"% var)"))
  
  
  pdf(paste0("PCA_plots/varying_fraction_top_peaks/PCA_CPM_merged_peaks_EV_p30_p42_LPS_UT.top",fr*100,"PercentMoreVarPeaks.pdf"),onefile = T)
  print(g1)
  print(g2)
  print(g3)
  dev.off()
  
}


sampleDists <- dist(scaled_counts)
sampleDists


sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
library(RColorBrewer)
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
dir.create("heatmaps")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 6)
pdf("heatmaps/pheatmap_sample_clustering_allpeaks.pdf")

print(p)
dev.off()



#log2
sampleDists <- dist(log2(scaled_counts+0.01))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
library(RColorBrewer)
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
dir.create("heatmaps")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 12)
pdf("heatmaps/pheatmap_sample_clustering_allpeaks.log2+0.01.pdf")

print(p)
dev.off()

manual_peaks=read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/macs2_pooled_reps/merged_peaks/merged_peaks_EV.p30.p42_LPS.UT.sorted.bed.withcolnames.annotated.csv")

peaks_in_2ormore=manual_peaks[manual_peaks$nsamples>1,]

cdf=countdata[rownames(countdata)%in%peaks_in_2ormore$peakID,]
cdflibs <- colSums(cdf)
scaled_counts <- apply(cdf,1, function(x)x/cdflibs*1000000)
pca <- prcomp(scaled_counts,scale. = T)
dtp <- as.data.frame(pca$x)
dtp$samples <- rownames(dtp)
dtp$vector=gsub("(.*)_(.*)_(.*)","\\1",dtp$samples)
dtp$vector[dtp$vector=="P42"] <- "p42"
dtp$LPS=gsub("(.*)_(.*)_(.*)","\\2",dtp$samples)
dtp$rep=gsub("(.*)_(.*)_(.*)","\\3",dtp$samples)

#percent variance:
percvar=round(pca$sdev^2/sum(pca$sdev^2)*100,2)

setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/")
dir.create("PCA_plots")

pdf("PCA_plots/PCA_CPM_merged_peaks_EV_p30_p42_LPS_UT.PeaksInAtleast2Samples.pdf",onefile = T)
ggplot(dtp, aes(x=PC1,y=PC2,col=vector)) + 
  geom_point() +theme_classic() + xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
ggplot(dtp, aes(x=PC1,y=PC2,col=LPS)) + 
  geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
ggplot(dtp, aes(x=PC1,y=PC2,col=rep)) + 
  geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
dev.off()


sampleDists <- dist(scaled_counts)
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
#colnames(sampleDistMatrix) <- NULL
library(RColorBrewer)
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
dir.create("heatmaps")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 6)
pdf("heatmaps/pheatmap_sample_clustering_PeaksInAtleast2Samples.pdf")

print(p)
dev.off()
#log2
sampleDists <- dist(log2(scaled_counts+0.01))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
library(RColorBrewer)
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
dir.create("heatmaps")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 6)
pdf("heatmaps/pheatmap_sample_clustering_PeaksInAtleast2Samples.log2+0.01.pdf")
print(p)
dev.off()


#select only peaks with a means score higer than...
plot(density(manual_peaks$meanscore))
summary(manual_peaks$meanscore)
mean(manual_peaks$meanscore[manual_peaks$nsamples==5])

#will take higher than the mean 21
pl=list(name="peaks_-logqvalabove21",
        filter=manual_peaks$peakID[manual_peaks$meanscore>21])


cdf=countdata[rownames(countdata)%in%pl$filter,]
cdflibs <- colSums(cdf)
scaled_counts <- apply(cdf,1, function(x)x/cdflibs*1000000)
pca <- prcomp(scaled_counts,scale. = T)
dtp <- as.data.frame(pca$x)
dtp$samples <- rownames(dtp)
dtp$vector=gsub("(.*)_(.*)_(.*)","\\1",dtp$samples)
dtp$vector[dtp$vector=="P42"] <- "p42"
dtp$LPS=gsub("(.*)_(.*)_(.*)","\\2",dtp$samples)
dtp$rep=gsub("(.*)_(.*)_(.*)","\\3",dtp$samples)

#percent variance:
percvar=round(pca$sdev^2/sum(pca$sdev^2)*100,2)

setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/")
dir.create("PCA_plots")

pdf(paste0("PCA_plots/PCA_CPM_merged_peaks_EV_p30_p42_LPS_UT.",pl$name,".pdf"),onefile = T)
ggplot(dtp, aes(x=PC1,y=PC2,col=vector)) + 
  geom_point() +theme_classic() + xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
ggplot(dtp, aes(x=PC1,y=PC2,col=LPS)) + 
  geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
ggplot(dtp, aes(x=PC1,y=PC2,col=rep)) + 
  geom_point(aes(shape=vector)) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))
dev.off()
pdf(paste0("PCA_plots/PCA_CPM_merged_peaks_EV_p30_p42_LPS_UT.",pl$name,"LPScolor.pdf"))
ggplot(dtp, aes(x=PC1,y=PC2,col=LPS)) + 
  geom_point(aes(shape=vector),size=4) +theme_classic()+ xlab(paste0("PC1(",percvar[1],"% var)"))+
  ylab(paste0("PC2(",percvar[2],"% var)"))+theme(text = element_text(size=20))
dev.off()

sampleDists <- dist(scaled_counts)
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
#colnames(sampleDistMatrix) <- NULL
library(RColorBrewer)
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
dir.create("heatmaps")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 6)
pdf(paste0("heatmaps/pheatmap_sample_clustering.",pl$name,".pdf"))

print(p)
dev.off()

sampleDists <- dist(log2(scaled_counts+0.01))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
library(RColorBrewer)
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
dir.create("heatmaps")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 12)
pdf(paste0("heatmaps/pheatmap_sample_clustering.",pl$name,".log2+0.01.pdf"))

print(p)
dev.off()
