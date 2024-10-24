options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
## ---------------------------
setwd("/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind/ER")
samples=read.csv("sample_sheet_ER.csv")
bsgmm39=readRDS("/scratch/llorenzi/BSgen.mm39.rds")

library(DiffBind)
samples
#remove bad samples:
bad_samples=c("p42_LPS_2_Apr22","p30_UT_1_Apr22")
samples=samples[!samples$SampleID%in%bad_samples,]

dir.create("no_bad_samples")
setwd("no_bad_samples")

Mdba=dba(sampleSheet=samples,minOverlap=1,bRemoveM=F)

Mdba

newfeat=paste0(Mdba$samples$Condition,Mdba$samples$Treatment)
Mdba$class[DBA_CONDITION,]=newfeat #this works


pdf(paste0("plotMdba.occupancy.ER.pdf"))
plot(Mdba)
dev.off()

pdf(paste0("plotMdba.PCA.occupancy.ER.pdf"))
dba.plotPCA(Mdba)
dev.off()

Mdba.GL=dba.blacklist(Mdba,blacklist = F,greylist = bsgmm39)
# > Mdba
# 9 Samples, 39442 sites in matrix:
#   ID Factor Condition Treatment Replicate Intervals
# 1  p30_UT_1_Nov21     ER     p30UT        UT         1     34129
# 2  p42_UT_1_Nov21     ER     p42UT        UT         1     23381
# 3 p30_LPS_1_Apr22     ER    p30LPS       LPS         1     16022
# 4  p42_UT_2_Apr22     ER     p42UT        UT         2     23381
# 5 p30_LPS_2_Apr22     ER    p30LPS       LPS         2     16022
# 6 p30_LPS_3_Apr22     ER    p30LPS       LPS         3     16022
# 7 p42_LPS_2_Apr22     ER    p42LPS       LPS         2      2879
# 8  p30_UT_1_Apr22     ER     p30UT        UT         1     34129
# 9 p42_LPS_1_Apr22     ER    p42LPS       LPS         1      2879
# > Mdba.GL=dba.blacklist(Mdba,blacklist = F,greylist = bsgmm39)
# Counting control reads for greylist...
# Building greylist: /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/27NOV_re_split_bam_files/INP_p30_ER.sorted.markedDups.bam
# coverage: 9072917 bp (0.33%)
# Building greylist: /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/27NOV_re_split_bam_files/INP_p42_ER.sorted.markedDups.bam
# coverage: 9144278 bp (0.34%)
# Building greylist: /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind/new_bam_files_may22/INP_p30_LPS_ER.sorted.markedDups.proper_pairs.minq2.bam
# coverage: 3159181 bp (0.12%)
# Building greylist: /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind/new_bam_files_may22/INP_p42_UT_ER.sorted.markedDups.proper_pairs.minq2.bam
# coverage: 3025037 bp (0.11%)
# Building greylist: /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind/new_bam_files_may22/INP_p42_LPS_ER.sorted.markedDups.proper_pairs.minq2.bam
# coverage: 2795661 bp (0.10%)
# Building greylist: /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind/new_bam_files_may22/INP_p30_UT_ER.sorted.markedDups.proper_pairs.minq2.bam
# coverage: 2901133 bp (0.11%)
# INP_p30_ER.nov21: 528 ranges, 9072917 bases
# INP_p42_ER.nov21: 535 ranges, 9144278 bases
# INP_p30_LPS_ER.apr22: 265 ranges, 3159181 bases
# INP_p42_UT_ER.apr22: 252 ranges, 3025037 bases
# INP_p42_LPS_ER.apr22: 240 ranges, 2795661 bases
# INP_p30_UT_ER.apr22: 242 ranges, 2901133 bases
# Master greylist: 625 ranges, 10250198 bases
# Removed: 3205 of 168844 intervals.
# Removed: 682 merged (of 39442) and 682 (of 39442) consensus.
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
saveRDS(Mdba.count,"Mdba.counts.ER.RDS")

pdf(paste0("CEBPa_affinity_based_heatmap.pdf"))
plot(Mdba.count)
dev.off()

pdf(paste0("CEBPa_affinity_based_PCA.pdf"))
dba.plotPCA(Mdba.count,  attributes=c(DBA_CONDITION), label=DBA_ID)
dev.off()
Mdba.count=readRDS("Mdba.counts.ER.RDS")

# #Establishing a contrast
# Mdba <- dba.contrast(Mdba, categories=DBA_FACTOR, minMembers = 2)
# 
#DBA_nov21=dba(Mdba.count,mask = c(1,2))
#normCounts <- dba.peakset(DBA_nov21, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
#head(normCounts)
# > head(normCounts)
# CHR   START     END p30_UT_1_Nov21 p42_UT_1_Nov21
# 1 chr1 3051797 3052197       62.05941       25.68196
# 2 chr1 3066404 3066804       59.71755       45.24916
# 3 chr1 3122831 3123231       87.81992       62.37047
# 4 chr1 3129802 3130202      539.79980      311.85235
# 5 chr1 3164881 3165281       53.86289       26.90491
# 6 chr1 3168151 3168551       58.54662       80.71473
#DBA_nov21 <- dba.contrast(DBA_nov21,design = '~ Condition')
# Computing results names...
# Error in DESeq2::estimateDispersionsGeneEst(res, quiet = TRUE, maxit = 10) : 
#   the number of samples and the number of model coefficients are equal,
# i.e., there are no replicates to estimate the dispersion.
# use an alternate design formula

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

ctsts=list("p30UT_vs_p42UT"=c(0,-1,0,0),   
           "p30LPS_vs_p42LPS"=c(0,0,1,-1),
           "p30LPS_vs_p30UT"=c(0,0,1,0),  
           "p42LPS_vs_p42UT"=c(0,-1,0,1) )

Dba=DBA
for (ct in ctsts) {
  Dba <- dba.contrast(Dba, contrast = ct)
  
}

Mdba.analysis=dba.analyze(Dba,method = DBA_ALL_METHODS)
saveRDS(Mdba.analysis,"Mdba.analysis.ER.RDS")

peaks=dba.peakset(Mdba.analysis,bRetrieve = T,
                  DataType = DBA_DATA_FRAME)

write.csv(peaks,"Consensus_peaks.DiffBind.ER.csv",row.names = F)

