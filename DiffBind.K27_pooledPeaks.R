## ---------------------------
##
##
## Purpose of script: Run DiffBind for CEBPa samples
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-04-25
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------
#module load conda/current
#conda activate atacseqqc
#R
options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
## ---------------------------
setwd("/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind")
samples=read.csv("sample_sheet_csuc.pooledReps.csv")
bsgmm39=readRDS("/scratch/llorenzi/BSgen.mm39.rds")

library(DiffBind)
Mdba=dba(sampleSheet=samples,minOverlap=1)

Mdba

newfeat=paste0(Mdba$samples$Condition,Mdba$samples$Treatment)
Mdba$class[DBA_CONDITION,]=newfeat #this works

dir.create("K27_pooledPeaks")
setwd("K27_pooledPeaks")

pdf(paste0("plotMdba.New_april22.pdf"))
plot(Mdba)
dev.off()

pdf(paste0("plotMdba.PCA.New_april22.pdf"))
dba.plotPCA(Mdba)
dev.off()

Mdba.GL=dba.blacklist(Mdba,blacklist = F,greylist = bsgmm39)
Mdba_greylist <- dba.blacklist(Mdba.GL, Retrieve=DBA_GREYLIST)
saveRDS(Mdba_greylist,"Mdba.greylist.RDS")

Mdba.GL=dba(Mdba.GL,mask=Mdba.GL$mask$All,minOverlap = 1)

#Affinity binding matrix
#The next step is to take the alignment files and compute 
#count information for each of the peaks/regions in the consensus set. 
#In this step, for each of the consensus regions DiffBind takes the number of 
#aligned reads in the ChIP sample and the input sample, to compute a normalized read 
#count for each sample at every potential binding site. The peaks in the consensus 
#peakset may be re-centered and trimmed based on calculating their summits 
#(point of greatest read overlap) in order to provide more standardized peak intervals.

# #bUseSummarizeOverlaps: logical indicating that <E2><80><98>summarizeOverlaps<E2><80><99>
#   should be used for counting instead of the built-in counting
# code.  This option is slower but uses the more standard
# counting function. If <E2><80><98>TRUE<E2><80><99>, all read files must be BAM
# (.bam extension), with associated index files (.bam.bai
#                                                extension).  The <E2><80><98>fragmentSize<E2><80><99> parameter must absent

#library(Rsamtools)
#Mdba$config$scanbamparam=ScanBamParam(mapqFilter = 2)
library(Rsamtools)
Mdba.GL$config$mapQCth <- 2
Mdba.GL$config$scanbamparam <- ScanBamParam(flag = scanBamFlag(isProperPair=T),mapqFilter=2)
Mdba.GL$config$fragmentSize=0 
Mdba.GL$config$singleEnd=F #actually, after reading the documentation
#I realise that setting these two things is not necessary when bUseSummarizeOverlaps=TRUE (default)
#because this counting mode "The fragmentSize parameter must absent." 
#AND DBA$config$singleEnd	logical indicating if reads are 
#single end; if NULL, status will be automatically detected.

Mdba.count=dba.count(Mdba.GL)
saveRDS(Mdba.count,"Mdba.counts.RDS")

pdf(paste0("CEBPa_affinity_based_heatmap.pdf"))
plot(Mdba.count)
dev.off()

pdf(paste0("CEBPa_affinity_based_PCA.pdf"))
dba.plotPCA(Mdba,  attributes=c(DBA_CONDITION,DBA_TREATMENT), label=DBA_ID)
dev.off()
Mdba.count=readRDS("Mdba.counts.RDS")

# #Establishing a contrast
# Mdba <- dba.contrast(Mdba, categories=DBA_FACTOR, minMembers = 2)
# 
DBA <- dba.contrast(Mdba.count,design = '~ Condition')
DBA$DESeq2$names
# [1] "Intercept"                 "Condition_EVUT_vs_EVLPS"  
# [3] "Condition_p30LPS_vs_EVLPS" "Condition_p30UT_vs_EVLPS" 
# [5] "Condition_p42LPS_vs_EVLPS" "Condition_p42UT_vs_EVLPS" 


#If I set the contrasts in the "easy" way, I get an error:
# ctsts <- list(c("Condition", "p30UT", "p42UT"),
#               c("Condition", "p30LPS", "p42LPS"),
#               c("Condition","p42LPS","EVLPS"),
#               c("Condition","p30LPS","EVLPS"),
#               c("Condition","p42UT","EVUT"),
#               c("Condition","p30UT","EVUT"),
#               c("Condition", "p30LPS", "p30UT"),
#               c("Condition", "p42LPS", "p42UT"),
#               c("Condition", "EVLPS", "EVUT"))
# Dba=DBA
# for (ct in ctsts) {
#   Dba <- dba.contrast(Dba, contrast = ct)
#   
# }
# Error in pv.contrastDesign(pv = pv, design = design, contrast = contrast,  : 
#                              Invalid contrast: no replicates in one group.


# "p30UT_vs_p42UT:0,0,0,1,0,-1"   
# "p30LPS_vs_p42LPS:0,0,1,0,-1,0"
# "p42LPS_vs_EVLPS:0,0,0,0,1,0"   
# "p30LPS_vs_EVLPS:0,0,1,0,0,0"  
# "p42UT_vs_EVUT:0,-1,0,0,0,1"    
# "p30UT_vs_EVUT:0,-1,0,1,0,0"   
# "p30LPS_vs_p30UT:0,0,1,-1,0,0"  
# "p42LPS_vs_p42UT:0,0,0,0,1,-1" 
# "EVLPS_vs_EVUT:0,-1,0,0,0,0" 

ctsts=list("p30UT_vs_p42UT"=c(0,0,0,1,0,-1),   
           "p30LPS_vs_p42LPS"=c(0,0,1,0,-1,0),
           "p42LPS_vs_EVLPS"=c(0,0,0,0,1,0),   
           "p30LPS_vs_EVLPS"=c(0,0,1,0,0,0),  
           "p42UT_vs_EVUT"=c(0,-1,0,0,0,1),    
           "p30UT_vs_EVUT"=c(0,-1,0,1,0,0),   
           "p30LPS_vs_p30UT"=c(0,0,1,-1,0,0),  
           "p42LPS_vs_p42UT"=c(0,0,0,0,1,-1), 
           "EVLPS_vs_EVUT"=c(0,-1,0,0,0,0) )

Dba=DBA
for (ct in ctsts) {
  Dba <- dba.contrast(Dba, contrast = ct)
  
}

Mdba.analysis=dba.analyze(Dba,method = DBA_ALL_METHODS)
saveRDS(Mdba.analysis,"Mdba.analysis.apr22.pooledReps.RDS")

peaks=dba.peakset(Mdba.analysis,bRetrieve = T,
                  DataType = DBA_DATA_FRAME)

write.csv(peaks,"onsensus_peaks.DiffBind.pooledReps.csv",row.names = F)

