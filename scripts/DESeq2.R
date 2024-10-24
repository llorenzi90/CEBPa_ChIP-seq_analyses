## ---------------------------
##
##
## Purpose of script:DESeq2 CEBEPa
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-02-28
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

#DESeq2 comparisons. We want to compare p30 vs p42 when treated with
#LPS and when not treated:
#based on the useful Deseq2 vignette; https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions
# dds$group <- factor(paste0(dds$genotype, dds$condition))
# design(dds) <- ~ group
# dds <- DESeq(dds)
# resultsNames(dds)
# results(dds, contrast=c("group", "IB", "IA"))
dtp$group <- factor(paste0(dtp$vector,dtp$LPS))
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = dtp,
                                 design = ~ group)
dds <- DESeq(ddsMat)
resultsNames(dds)
res1=results(dds, contrast=c("group", "p30UT", "p42UT"))
summary(res1)
res2=results(dds, contrast=c("group", "p30LPS", "p42LPS"))
summary(res2)
summary(results(dds,contrast = c("group","p30LPS","EVLPS")))
summary(results(dds,contrast = c("group","p42LPS","EVLPS")))
summary(results(dds,contrast = c("group","p30UT","EVUT")))
summary(results(dds,contrast = c("group","p42UT","EVUT")))
# 
# ctsts <- list(c("group", "p30UT", "p42UT"),
#               c("group", "p30LPS", "p42LPS"),
#               c("group","p42LPS","EVLPS"),
#               c("group","p30LPS","EVLPS"),
#               c("group","p42UT","EVUT"),
#               c("group","p30UT","EVUT"))

#17 marzo: p30LPSvsp30UT, p42LPSvsp42UT y EVLPSvsEV UT
ctsts <- list(c("group", "p30LPS", "p30UT"),
              c("group", "p42LPS", "p42UT"),
              c("group", "EVLPS", "EVUT"))
dir.create("DESeq2_results")
setwd("DESeq2_results/")
for (cts in ctsts) {
  res=results(dds, contrast=cts)
  nam=paste0(cts[2],"vs",cts[3])
  print(summary(res))
  write(capture.output(summary(res)),paste0(nam,"_summary_res.txt"))
  
  res <- as.data.frame(res)
  res <- res[order(res$padj),]
  res <- res[!is.na(res$padj),]
  write.csv(res,paste0(nam,"_DESeq2_res.noNAs.csv"))
  
}
write.csv(dtp[,13:17],"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/H3K27ac_sample_metadata.txt",row.names = F)