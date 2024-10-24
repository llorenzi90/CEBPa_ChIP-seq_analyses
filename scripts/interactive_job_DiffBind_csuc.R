Mdba <- readRDS("DBA.keptlength.RDS")
Mdba2=Mdba
newfeat=paste0(Mdba$samples$Condition,Mdba$samples$Treatment)
Mdba$samples$Condition=newfeat #this does not work
Mdba$class[DBA_CONDITION,]=newfeat #this does
DBA <- dba.contrast(Mdba, 
                    design = '~ Condition')
#DBA2 <- dba.contrast(Mdba2, design = '~ Condition * Treatment',reorderMeta=list(Condition=c("p42","p30","EV"),
                                                                              
ctsts <- list(c("Condition", "p30UT", "p42UT"),
                            c("Condition", "p30LPS", "p42LPS"),
                            c("Condition","p42LPS","EVLPS"),
                            c("Condition","p30LPS","EVLPS"),
                            c("Condition","p42UT","EVUT"),
                            c("Condition","p30UT","EVUT"),
                            c("Condition", "p30LPS", "p30UT"),
                            c("Condition", "p42LPS", "p42UT"),
                            c("Condition", "EVLPS", "EVUT"))
Dba=DBA
for (ct in ctsts) {
  Dba <- dba.contrast(Dba, contrast = ct)
  
}

ctsts2 <- list(c(0,0,0,1,0,-1), c(0,0,1,0,-1,0),c(0,0,0,0,1,0))
#c("Condition", "p30UT", "p42UT")== 0,0,0,1,0,-1
#c("Condition", "p30LPS", "p42LPS")==c(0,0,1,0,-1,0)
#c("Condition","p42LPS","EVLPS") == c(0,0,0,0,1,0)

Dbatest=DBA
for (ct in ctsts2) {
  Dbatest <- dba.contrast(Dbatest, contrast = ct)
  
}

Dba=dba.analyze(Dba)
Dbatest=dba.analyze(Dbatest)

res1 <- lapply(Dba$contrasts, function(x){
  re=x$DESeq2$de
  re$diff=ifelse(re$padj<=0.05,ifelse(re$fold>0,"UP","DOWN"),"NO")
  return(re)
})

res2 <- lapply(Dbatest$contrasts, function(x){
  re=x$DESeq2$de
  re$diff=ifelse(re$padj<=0.05,ifelse(re$fold>0,"UP","DOWN"),"NO")
  return(re)
})


lapply(res1, function(x)table(x$diff))
lapply(res2, function(x)table(x$diff))

lapply(res1, function(x)table(x$diff))
# [[1]]
# 
# DOWN     NO     UP 
# 5926 159038   3645 
# 
# [[2]]
# 
# DOWN     NO     UP 
# 4484 156521   7604 

lapply(res2, function(x)table(x$diff))
# [[1]]
# 
# DOWN     NO     UP 
# 5926 159038   3645 
# 
# [[2]]
# 
# DOWN     NO     UP 
# 4484 156521   7604 
# 
# [[3]]
# 
# DOWN     NO     UP 
# 5402 157199   6008 


dba.show(Dba, bContrast=T)
# Factor  Group Samples Group2 Samples2 DB.DESeq2
# 1 Condition  p30UT       3  p42UT        2      9571
# 2 Condition p30LPS       3 p42LPS        2     12088
dba.show(Dbatest, bContrast=T)
# Factor        Group Intercept Condition_EVUT_vs_EVLPS
# 1 Coefficient 0,0,0,1,0,-1         0                       0
# 2 Coefficient 0,0,1,0,-1,0         0                       0
# 3 Coefficient  0,0,0,0,1,0         0                       0
# Condition_p30LPS_vs_EVLPS Condition_p30UT_vs_EVLPS Condition_p42LPS_vs_EVLPS
# 1                         0                        1                         0
# 2                         1                        0                        -1
# 3                         0                        0                         1
# Condition_p42UT_vs_EVLPS DB.DESeq2
# 1                       -1      9571
# 2                        0     12088
# 3                        0     11410

Dba_res=lapply(1:length(Dba$contrasts),function(x)
  dba.report(Dba,contrast = x,DataType = 'DBA_DATA_FRAME'))

c("Condition","p30LPS","EVLPS"),
c("Condition","p42UT","EVUT"),
c("Condition","p30UT","EVUT"),
c("Condition", "p30LPS", "p30UT"),
c("Condition", "p42LPS", "p42UT"),
c("Condition", "EVLPS", "EVUT"))

saveRDS(Dbatest,"Dba_first_3comps.RDS")
remaining_ctsts=list(c(0,0,1,0,0,0),
                     c(0,-1,0,0,0,1),
                     c(0,-1,0,1,0,0),
                     c(0,0,1,-1,0,0),
                     c(0,0,0,0,1,-1),
                     c(0,-1,0,0,0,0))

newdba=DBA

for (ct in remaining_ctsts) {
  newdba=dba.contrast(newdba,contrast = ct)
}

remaining_ctsts=list(c(0,0,1,0,0,0),
                     c(0,-1,0,0,0,1),
                     c(0,-1,0,1,0,0),
                     c(0,0,1,-1,0,0),
                     c(0,0,0,0,1,-1),
                     c(0,-1,0,0,0,0))


# 
# [llorenzi@ijc20844 ~]$ sudo mount shares/INVESTIGACIO/
#   [llorenzi@ijc20844 ~]$ csuc
# Last login: Fri Mar 18 23:01:54 2022 from 192.168.79.1
# +-----------------------------------------------------------------------+
#   |                      Welcome to CSUC HPC Service                      |
#   +-----------------------------------------------------------------------+
#   |            Use 'module av' to see available apps and tools            |
#   |                                                                       |
#   |               Runtimes have to be specified for all jobs              |
#   |                                                                       |
#   |   Check our portal for documentation & support: http://hpc.csuc.cat   |
#   +-----------------------------------------------------------------------+
#   +-----------+-------------+-----------------+--------------+------------+
#   |  MACHINE  | TOTAL SLOTS | ALLOCATED SLOTS | QUEUED SLOTS | OCCUPATION |
#   +-----------+-------------+-----------------+--------------+------------+
#   | std nodes |        2016 |            1760 |          820 |       87 % |
#   | fat nodes |         288 |             283 |            0 |       98 % |
#   | mem nodes |         384 |             299 |          160 |       77 % |
#   | gpu nodes |         192 |             144 |          336 |       75 % |
#   | knl nodes |        1088 |               0 |            0 |        0 % |
#   | res nodes |         864 |             523 |            0 |       60 % |
#   +-----------+-------------+-----------------+--------------+------------+
#   +-----------------------------------------------------------------------+
#   |                        CSUC HPC Service: consum                       |
#   |                     UCs consumed by user llorenzi                     |
#   |                     from 2022-01-01 to 2023-01-01                     |
#   +-----------------------------------------------------------------------+
#   |  Account                  llorenzi            Group         Assigned  |
#   +-----------------------------------------------------------------------+
#   |  ijcres_knl                    0.0              0.0           100000  |
#   |  ijcres_low                    0.0              0.0                0  |
#   |  ijcres_normal              4480.7          11144.6            50000  |
#   +-----------------------------------------------------------------------+
#   See "consum --help" for more details.
# llorenzi@login2:/home/llorenzi>salloc -t 120 --ntasks=10
# salloc: Pending job allocation 1962228
# salloc: job 1962228 queued and waiting for resources
# salloc: error: Lookup failed: Unknown host
# salloc: job 1962228 has been allocated resources
# salloc: Granted job allocation 1962228
# Your SCRATCH directory is located at pirineusgpu1 /tmp/llorenzi/1962228
# and will be automatically copied to /scratch/llorenzi/tmp/1962228
# llorenzi@pirineusgpu1:/home/llorenzi>cd /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind/
#   llorenzi@pirineusgpu1:/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind>module load conda/current 
# WARNING: Note that a configuration file has been sourced and won't be reversed by module unload or module switch
# llorenzi@pirineusgpu1:/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind>conda activate atacseqqc
# (atacseqqc) llorenzi@pirineusgpu1:/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind>R
# 
# R version 4.1.1 (2021-08-10) -- "Kick Things"
# Copyright (C) 2021 The R Foundation for Statistical Computing
# Platform: x86_64-conda-linux-gnu (64-bit)
# 
# R es un software libre y viene sin GARANTIA ALGUNA.
# Usted puede redistribuirlo bajo ciertas circunstancias.
# Escriba 'license()' o 'licence()' para detalles de distribucion.
# 
# R es un proyecto colaborativo con muchos contribuyentes.
# Escriba 'contributors()' para obtener más información y
# 'citation()' para saber cómo citar R o paquetes de R en publicaciones.
# 
# Escriba 'demo()' para demostraciones, 'help()' para el sistema on-line de ayuda,
# o 'help.start()' para abrir el sistema de ayuda HTML con su navegador.
# Escriba 'q()' para salir de R.
# 
# > library(DiffBind)
# Loading required package: GenomicRanges
# Loading required package: stats4
# Loading required package: BiocGenerics
# Loading required package: parallel
# 
# Attaching package: ‘BiocGenerics’
# 
# The following objects are masked from ‘package:parallel’:
# 
#     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
#     clusterExport, clusterMap, parApply, parCapply, parLapply,
#     parLapplyLB, parRapply, parSapply, parSapplyLB
# 
# The following objects are masked from ‘package:stats’:
# 
#     IQR, mad, sd, var, xtabs
# 
# The following objects are masked from ‘package:base’:
# 
#     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
#     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
#     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
#     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
#     union, unique, unsplit, which.max, which.min
# 
# Loading required package: S4Vectors
# 
# Attaching package: ‘S4Vectors’
# 
# The following objects are masked from ‘package:base’:
# 
#     expand.grid, I, unname
# 
# Loading required package: IRanges
# Loading required package: GenomeInfoDb
# Loading required package: SummarizedExperiment
# Loading required package: MatrixGenerics
# Loading required package: matrixStats
# 
# Attaching package: ‘MatrixGenerics’
# 
# The following objects are masked from ‘package:matrixStats’:
# 
#     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#     colWeightedMeans, colWeightedMedians, colWeightedSds,
#     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#     rowWeightedSds, rowWeightedVars
# 
# Loading required package: Biobase
# Welcome to Bioconductor
# 
#     Vignettes contain introductory material; view with
#     'browseVignettes()'. To cite Bioconductor, see
#     'citation("Biobase")', and for packages 'citation("pkgname")'.
# 
# 
# Attaching package: ‘Biobase’
# 
# The following object is masked from ‘package:MatrixGenerics’:
# 
#     rowMedians
# 
# The following objects are masked from ‘package:matrixStats’:
# 
#     anyMissing, rowMedians
# 
# 
#  >>> DiffBind 3.2.
# > Mdba <- readRDS("DBA.keptlength.RDS")
# 
# > 
# > Mdba2=Mdba
# > newfeat=paste0(Mdba$samples$Condition,Mdba$samples$Treatment)
# Mdba$samples$Condition=newfeat #this does not work
# Mdba$class[DBA_CONDITION,]=newfeat
# > DBA <- dba.contrast(Mdba, 
#                     design = '~ Condition')
# Computing results names...
# > DBA2=dba.contrast(Mdba2, design='~ Condition * Treatment')
# Computing results names...
# > ctsts <- list(c("Condition", "p30UT", "p42UT"),
#                             c("Condition", "p30LPS", "p42LPS"),
#                             c("Condition","p42LPS","EVLPS"),
#                             c("Condition","p30LPS","EVLPS"),
#                             c("Condition","p42UT","EVUT"),
#                             c("Condition","p30UT","EVUT"),
#                             c("Condition", "p30LPS", "p30UT"),
#                             c("Condition", "p42LPS", "p42UT"),
#                             c("Condition", "EVLPS", "EVUT"))
# > Dba=DBA
# for (ct in ctsts) {
#   Dba <- dba.contrast(Dba, contrast = ct)
#   
# }
# Error in pv.contrastDesign(pv = pv, design = design, contrast = contrast,  : 
#   Invalid contrast: no replicates in one group.
# > ct
# [1] "Condition" "p42LPS"    "EVLPS"    
# > test=dba.contrast(DBA2,contrast=list(c(0,0,0,0,1,0)))
# Error: Invalid contrast. Names in list must match one of design column names.
# Además: Warning message:
# In if (!contrast[[1]] %in% names) { :
#   la condición tiene longitud > 1 y sólo el primer elemento será usado
# > test=dba.contrast(DBA,contrast=list(c(0,0,0,0,1,0)))
# Error: Invalid contrast. Names in list must match one of design column names.
# Además: Warning message:
# In if (!contrast[[1]] %in% names) { :
#   la condición tiene longitud > 1 y sólo el primer elemento será usado
# > test=dba.contrast(DBA,contrast=c(0,0,0,0,1,0))
# > test$contrasts$
# test$contrasts$
# > test$contrasts
# test$contrasts
# > test$contrasts[[1]]$contrast
# [1] 0 0 0 0 1 0
# > test$contrasts[[1]]$contrast
# test$contrasts[[1]]$contrastType  test$contrasts[[1]]$contrast
# > test$contrasts[[1]]$contrastType
# [1] "bycolumn"
# > Dba$
# Dba$peaks        Dba$called       Dba$attributes   Dba$filterFun    Dba$DESeq2
# Dba$class        Dba$score        Dba$minOverlap   Dba$minCount     Dba$contrasts
# Dba$chrmap       Dba$binding      Dba$masks        Dba$summits      
# Dba$config       Dba$merged       Dba$SN           Dba$meta         
# Dba$samples      Dba$totalMerged  Dba$maxFilter    Dba$design       
# > Dba$con
# Dba$config     Dba$contrasts  
# > length(Dba$contrasts)
# [1] 2
# > Dba=dba.analyze(Dba)
# Applying Blacklist/Greylists...
# No genome detected.
# Normalize DESeq2 with defaults...
# Analyzing...
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# > Dba$con
# Dba$config     Dba$contrasts  
# > Dba$contrasts[[2]]$
# Dba$contrasts[[2]]$contrastType  Dba$contrasts[[2]]$group1
# Dba$contrasts[[2]]$contrast      Dba$contrasts[[2]]$group2
# Dba$contrasts[[2]]$name1         Dba$contrasts[[2]]$DESeq2
# Dba$contrasts[[2]]$name2         
# > head(Dba$contrasts[[2]]$DESeq2$de)
#            id         pval         padj      fold
# 87525   87525 1.794090e-62 2.966330e-57 -4.297844
# 54174   54174 2.312965e-56 1.912116e-51  3.875301
# 124964 124964 1.618434e-55 8.919676e-51 -3.428487
# 2056     2056 1.357360e-54 5.610613e-50  3.999920
# 81191   81191 1.950661e-46 6.450408e-42 -4.250808
# 157487 157487 7.326973e-46 2.019057e-41 -3.031725
# > Dba$contrasts[[1]]$group1
#  EV_LPS_R1   EV_UT_R1 p30_LPS_R1 p30_LPS_R2 p30_LPS_R3  p30_UT_R1  p30_UT_R2 
#      FALSE      FALSE      FALSE      FALSE      FALSE       TRUE       TRUE 
#  p30_UT_R3 p42_LPS_R1 p42_LPS_R2  p42_UT_R1  p42_UT_R2 
#       TRUE      FALSE      FALSE      FALSE      FALSE 
# > salloc: Job 1962228 has exceeded its time limit and its allocation has been revoked.
# srun: Job step aborted: Waiting up to 32 seconds for job step to finish.
# srun: error: pirineusgpu1: task 0: Killed
# llorenzi@login2:/home/llorenzi>salloc -t 120 --ntasks=10
# salloc: Pending job allocation 1962275
# salloc: job 1962275 queued and waiting for resources
#  ^Csalloc: Job allocation 1962275 has been revoked.
# salloc: Job aborted due to signal
# llorenzi@login2:/home/llorenzi>salloc -p mem -t 120 --ntasks=10
# salloc: Granted job allocation 1962276
# Your SCRATCH directory is located at canigo1 /tmp/llorenzi/1962276
# and will be automatically copied to /scratch/llorenzi/tmp/1962276
# llorenzi@canigo1:/home/llorenzi>cd /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/
# llorenzi@canigo1:/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples>cd DiffBind/
# llorenzi@canigo1:/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind>module load conda/current
# WARNING: Note that a configuration file has been sourced and won't be reversed by module unload or module switch
# llorenzi@canigo1:/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind>conda activate atacseqqc
# (atacseqqc) llorenzi@canigo1:/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind>R
# 
# R version 4.1.1 (2021-08-10) -- "Kick Things"
# Copyright (C) 2021 The R Foundation for Statistical Computing
# Platform: x86_64-conda-linux-gnu (64-bit)
# 
# R es un software libre y viene sin GARANTIA ALGUNA.
# Usted puede redistribuirlo bajo ciertas circunstancias.
# Escriba 'license()' o 'licence()' para detalles de distribucion.
# 
# R es un proyecto colaborativo con muchos contribuyentes.
# Escriba 'contributors()' para obtener más información y
# 'citation()' para saber cómo citar R o paquetes de R en publicaciones.
# 
# Escriba 'demo()' para demostraciones, 'help()' para el sistema on-line de ayuda,
# o 'help.start()' para abrir el sistema de ayuda HTML con su navegador.
# Escriba 'q()' para salir de R.
# 
# > library(DiffBind)
# Loading required package: GenomicRanges
# Loading required package: stats4
# Loading required package: BiocGenerics
# Loading required package: parallel
# 
# Attaching package: ‘BiocGenerics’
# 
# The following objects are masked from ‘package:parallel’:
#   
#   clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
# clusterExport, clusterMap, parApply, parCapply, parLapply,
# parLapplyLB, parRapply, parSapply, parSapplyLB
# 
# The following objects are masked from ‘package:stats’:
#   
#   IQR, mad, sd, var, xtabs
# 
# The following objects are masked from ‘package:base’:
#   
#   anyDuplicated, append, as.data.frame, basename, cbind, colnames,
# dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
# grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
# order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
# rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
# union, unique, unsplit, which.max, which.min
# 
# Loading required package: S4Vectors
# 
# Attaching package: ‘S4Vectors’
# 
# The following objects are masked from ‘package:base’:
#   
#   expand.grid, I, unname
# 
# Loading required package: IRanges
# Loading required package: GenomeInfoDb
# Loading required package: SummarizedExperiment
# Loading required package: MatrixGenerics
# Loading required package: matrixStats
# 
# Attaching package: ‘MatrixGenerics’
# 
# The following objects are masked from ‘package:matrixStats’:
#   
#   colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
# colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
# colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
# colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
# colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
# colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
# colWeightedMeans, colWeightedMedians, colWeightedSds,
# colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
# rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
# rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
# rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
# rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
# rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
# rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
# rowWeightedSds, rowWeightedVars
# 
# Loading required package: Biobase
# Welcome to Bioconductor
# 
# Vignettes contain introductory material; view with
# 'browseVignettes()'. To cite Bioconductor, see
# 'citation("Biobase")', and for packages 'citation("pkgname")'.
# 
# 
# Attaching package: ‘Biobase’
# 
# The following object is masked from ‘package:MatrixGenerics’:
#   
#   rowMedians
# 
# The following objects are masked from ‘package:matrixStats’:
#   
#   anyMissing, rowMedians
# 
# 
# >>> DiffBind 3.2.
# > newfeat=paste0(Mdba$samples$Condition,Mdba$samples$Treatment)
# Mdba$samples$Condition=newfeat #this does not work
# Mdba$class[DBA_CONDITION,]=newfeat #this does
# 
# > Mdba <- readRDS("DBA.keptlength.RDS")
# Mdba2=Mdba
# newfeat=paste0(Mdba$samples$Condition,Mdba$samples$Treatment)
# Mdba$samples$Condition=newfeat #this does not work
# Mdba$class[DBA_CONDITION,]=newfeat #this does
# DBA <- dba.contrast(Mdba, 
#                     design = '~ Condition')
# DBA2 <- dba.contrast(Mdba2, design = '~ Condition * Treatment',reorderMeta=list(Condition=c("p42","p30","EV"),
#                                                                                 
#                                                                                 ctsts <- list(c("Condition", "p30UT", "p42UT"),
#                                                                                               c("Condition", "p30LPS", "p42LPS"),
#                                                                                               c("Condition","p42LPS","EVLPS"),
#                                                                                               c("Condition","p30LPS","EVLPS"),
#                                                                                               c("Condition","p42UT","EVUT"),
#                                                                                               c("Condition","p30UT","EVUT"),
#                                                                                               c("Condition", "p30LPS", "p30UT"),
#                                                                                               c("Condition", "p42LPS", "p42UT"),
#                                                                                               c("Condition", "EVLPS", "EVUT"))
#                                                                                 Dba=DBA
#                                                                                 for (ct in ctsts) {
#                                                                                   Dba <- dba.contrast(Dba, contrast = ct)
#                                                                                   
#                                                                                 }
#                                                                                 
#                                                                                 Computing results names...
#                                                                                 Error: unexpected symbol in:
#                                                                                   "                            c("Condition", "EVLPS", "EVUT"))
# Dba"
#                                                                                 > 
#                                                                                   > Mdba <- readRDS("DBA.keptlength.RDS")
#                                                                                 Mdba2=Mdba
#                                                                                 newfeat=paste0(Mdba$samples$Condition,Mdba$samples$Treatment)
#                                                                                 Mdba$samples$Condition=newfeat #this does not work
#                                                                                 Mdba$class[DBA_CONDITION,]=newfeat #this does
#                                                                                 DBA <- dba.contrast(Mdba, 
#                                                                                                     design = '~ Condition')
#                                                                                 DBA2 <- dba.contrast(Mdba2, design = '~ Condition * Treatment',reorderMeta=list(Condition=c("p42","p30","EV"),
#                                                                                                                                                                 
#                                                                                                                                                                 ctsts <- list(c("Condition", "p30UT", "p42UT"),
#                                                                                                                                                                               c("Condition", "p30LPS", "p42LPS"),
#                                                                                                                                                                               c("Condition","p42LPS","EVLPS"),
#                                                                                                                                                                               c("Condition","p30LPS","EVLPS"),
#                                                                                                                                                                               c("Condition","p42UT","EVUT"),
#                                                                                                                                                                               c("Condition","p30UT","EVUT"),
#                                                                                                                                                                               c("Condition", "p30LPS", "p30UT"),
#                                                                                                                                                                               c("Condition", "p42LPS", "p42UT"),
#                                                                                                                                                                               c("Condition", "EVLPS", "EVUT"))
#                                                                                                                                                                 Dba=DBA
#                                                                                                                                                                 for (ct in ctsts) {
#                                                                                                                                                                   Dba <- dba.contrast(Dba, contrast = ct)
#                                                                                                                                                                   
#                                                                                                                                                                 }
#                                                                                                                                                                 Computing results names...
#                                                                                                                                                                 Error: unexpected symbol in:
#                                                                                                                                                                   "                            c("Condition", "EVLPS", "EVUT"))
# Dba"
#                                                                                                                                                                 > ctsts
#                                                                                                                                                                 Error: objeto 'ctsts' no encontrado
#                                                                                                                                                                 > ctsts <- list(c("Condition", "p30UT", "p42UT"),
#                                                                                                                                                                                 c("Condition", "p30LPS", "p42LPS"),
#                                                                                                                                                                                 c("Condition","p42LPS","EVLPS"),
#                                                                                                                                                                                 c("Condition","p30LPS","EVLPS"),
#                                                                                                                                                                                 c("Condition","p42UT","EVUT"),
#                                                                                                                                                                                 c("Condition","p30UT","EVUT"),
#                                                                                                                                                                                 c("Condition", "p30LPS", "p30UT"),
#                                                                                                                                                                                 c("Condition", "p42LPS", "p42UT"),
#                                                                                                                                                                                 c("Condition", "EVLPS", "EVUT"))
#                                                                                                                                                                 > ctsts
#                                                                                                                                                                 [[1]]
#                                                                                                                                                                 [1] "Condition" "p30UT"     "p42UT"    
#                                                                                                                                                                 
#                                                                                                                                                                 [[2]]
#                                                                                                                                                                 [1] "Condition" "p30LPS"    "p42LPS"   
#                                                                                                                                                                 
#                                                                                                                                                                 [[3]]
#                                                                                                                                                                 [1] "Condition" "p42LPS"    "EVLPS"    
#                                                                                                                                                                 
#                                                                                                                                                                 [[4]]
#                                                                                                                                                                 [1] "Condition" "p30LPS"    "EVLPS"    
#                                                                                                                                                                 
#                                                                                                                                                                 [[5]]
#                                                                                                                                                                 [1] "Condition" "p42UT"     "EVUT"     
#                                                                                                                                                                 
#                                                                                                                                                                 [[6]]
#                                                                                                                                                                 [1] "Condition" "p30UT"     "EVUT"     
#                                                                                                                                                                 
#                                                                                                                                                                 [[7]]
#                                                                                                                                                                 [1] "Condition" "p30LPS"    "p30UT"    
#                                                                                                                                                                 
#                                                                                                                                                                 [[8]]
#                                                                                                                                                                 [1] "Condition" "p42LPS"    "p42UT"    
#                                                                                                                                                                 
#                                                                                                                                                                 [[9]]
#                                                                                                                                                                 [1] "Condition" "EVLPS"     "EVUT"     
#                                                                                                                                                                 
#                                                                                                                                                                 > Mdba <- readRDS("DBA.keptlength.RDS")
#                                                                                                                                                                 Mdba2=Mdba
#                                                                                                                                                                 newfeat=paste0(Mdba$samples$Condition,Mdba$samples$Treatment)
#                                                                                                                                                                 Mdba$samples$Condition=newfeat #this does not work
#                                                                                                                                                                 Mdba$class[DBA_CONDITION,]=newfeat #this does
#                                                                                                                                                                 DBA <- dba.contrast(Mdba, 
#                                                                                                                                                                                     design = '~ Condition')
#                                                                                                                                                                 DBA2 <- dba.contrast(Mdba2, design = '~ Condition * Treatment',reorderMeta=list(Condition=c("p42","p30","EV"),
#                                                                                                                                                                                                                                                 
#                                                                                                                                                                                                                                                 ctsts <- list(c("Condition", "p30UT", "p42UT"),
#                                                                                                                                                                                                                                                               c("Condition", "p30LPS", "p42LPS"),
#                                                                                                                                                                                                                                                               c("Condition","p42LPS","EVLPS"),
#                                                                                                                                                                                                                                                               c("Condition","p30LPS","EVLPS"),
#                                                                                                                                                                                                                                                               c("Condition","p42UT","EVUT"),
#                                                                                                                                                                                                                                                               c("Condition","p30UT","EVUT"),
#                                                                                                                                                                                                                                                               c("Condition", "p30LPS", "p30UT"),
#                                                                                                                                                                                                                                                               c("Condition", "p42LPS", "p42UT"),
#                                                                                                                                                                                                                                                               c("Condition", "EVLPS", "EVUT"))
#                                                                                                                                                                                                                                                 Dba=DBA
#                                                                                                                                                                                                                                                 for (ct in ctsts) {
#                                                                                                                                                                                                                                                   Dba <- dba.contrast(Dba, contrast = ct)
#                                                                                                                                                                                                                                                   
#                                                                                                                                                                                                                                                 }
#                                                                                                                                                                                                                                                 > DBA$
#                                                                                                                                                                                                                                                   DBA$peaks        DBA$called       DBA$attributes   DBA$filterFun    DBA$DESeq2
#                                                                                                                                                                                                                                                 DBA$class        DBA$score        DBA$minOverlap   DBA$minCount     
#                                                                                                                                                                                                                                                 DBA$chrmap       DBA$binding      DBA$masks        DBA$summits      
#                                                                                                                                                                                                                                                 DBA$config       DBA$merged       DBA$SN           DBA$meta         
#                                                                                                                                                                                                                                                 DBA$samples      DBA$totalMerged  DBA$maxFilter    DBA$design       
#                                                                                                                                                                                                                                                 > Dba=DBA
#                                                                                                                                                                                                                                                 for (ct in ctsts) {
#                                                                                                                                                                                                                                                   Dba <- dba.contrast(Dba, contrast = ct)
#                                                                                                                                                                                                                                                   
#                                                                                                                                                                                                                                                 }
#                                                                                                                                                                                                                                                 Error in pv.contrastDesign(pv = pv, design = design, contrast = contrast,  : 
#                                                                                                                                                                                                                                                                              Invalid contrast: no replicates in one group.
#                                                                                                                                                                                                                                                                            > DBA$
#                                                                                                                                                                                                                                                                              DBA$peaks        DBA$called       DBA$attributes   DBA$filterFun    DBA$DESeq2
#                                                                                                                                                                                                                                                                            DBA$class        DBA$score        DBA$minOverlap   DBA$minCount     
#                                                                                                                                                                                                                                                                            DBA$chrmap       DBA$binding      DBA$masks        DBA$summits      
#                                                                                                                                                                                                                                                                            DBA$config       DBA$merged       DBA$SN           DBA$meta         
#                                                                                                                                                                                                                                                                            DBA$samples      DBA$totalMerged  DBA$maxFilter    DBA$design       
#                                                                                                                                                                                                                                                                            > DBA$DESeq2
#                                                                                                                                                                                                                                                                            $names
#                                                                                                                                                                                                                                                                            [1] "Intercept"                 "Condition_EVUT_vs_EVLPS"  
#                                                                                                                                                                                                                                                                            [3] "Condition_p30LPS_vs_EVLPS" "Condition_p30UT_vs_EVLPS" 
#                                                                                                                                                                                                                                                                            [5] "Condition_p42LPS_vs_EVLPS" "Condition_p42UT_vs_EVLPS" 
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            > ctsts2 <- list(c(0,0,0,1,0,-1), c(0,0,1,0,-1,0))
#                                                                                                                                                                                                                                                                            > Dbatest=DBA
#                                                                                                                                                                                                                                                                            > ctsts2 <- list(c(0,0,0,1,0,-1), c(0,0,1,0,-1,0),c(0,0,0,0,1,0))
#                                                                                                                                                                                                                                                                            > Dbatest=DBA
#                                                                                                                                                                                                                                                                            for (ct in ctsts2) {
#                                                                                                                                                                                                                                                                              Dbatest <- dba.contrast(Dbatest, contrast = ct)
#                                                                                                                                                                                                                                                                              
#                                                                                                                                                                                                                                                                            }
#                                                                                                                                                                                                                                                                            > length(Dba$contrasts)
#                                                                                                                                                                                                                                                                            [1] 2
#                                                                                                                                                                                                                                                                            > length(Dbatest$contrasts)
#                                                                                                                                                                                                                                                                            [1] 3
#                                                                                                                                                                                                                                                                            > 
#                                                                                                                                                                                                                                                                              > Dba=dba.analyze(Dba)
#                                                                                                                                                                                                                                                                            Applying Blacklist/Greylists...
#                                                                                                                                                                                                                                                                            No genome detected.
#                                                                                                                                                                                                                                                                            Normalize DESeq2 with defaults...
#                                                                                                                                                                                                                                                                            Analyzing...
#                                                                                                                                                                                                                                                                            gene-wise dispersion estimates
#                                                                                                                                                                                                                                                                            mean-dispersion relationship
#                                                                                                                                                                                                                                                                            final dispersion estimates
#                                                                                                                                                                                                                                                                            > Dbatest=dba.analyze(Dbatest)
#                                                                                                                                                                                                                                                                            Applying Blacklist/Greylists...
#                                                                                                                                                                                                                                                                            No genome detected.
#                                                                                                                                                                                                                                                                            Normalize DESeq2 with defaults...
#                                                                                                                                                                                                                                                                            Analyzing...
#                                                                                                                                                                                                                                                                            gene-wise dispersion estimates
#                                                                                                                                                                                                                                                                            mean-dispersion relationship
#                                                                                                                                                                                                                                                                            final dispersion estimates
#                                                                                                                                                                                                                                                                            > Dba$
#                                                                                                                                                                                                                                                                              Dba$peaks        Dba$called       Dba$minOverlap   Dba$minCount     Dba$contrasts
#                                                                                                                                                                                                                                                                            Dba$class        Dba$binding      Dba$masks        Dba$summits      Dba$norm
#                                                                                                                                                                                                                                                                            Dba$chrmap       Dba$merged       Dba$SN           Dba$meta         Dba$score
#                                                                                                                                                                                                                                                                            Dba$config       Dba$totalMerged  Dba$maxFilter    Dba$design       
#                                                                                                                                                                                                                                                                            Dba$samples      Dba$attributes   Dba$filterFun    Dba$DESeq2       
#                                                                                                                                                                                                                                                                            > Dba$contrasts[[1]]$
#                                                                                                                                                                                                                                                                              Dba$contrasts[[1]]$contrastType  Dba$contrasts[[1]]$group1
#                                                                                                                                                                                                                                                                            Dba$contrasts[[1]]$contrast      Dba$contrasts[[1]]$group2
#                                                                                                                                                                                                                                                                            Dba$contrasts[[1]]$name1         Dba$contrasts[[1]]$DESeq2
#                                                                                                                                                                                                                                                                            Dba$contrasts[[1]]$name2         
#                                                                                                                                                                                                                                                                            > res1 <- lapply(Dba$contrasts, function(x){
#                                                                                                                                                                                                                                                                              re=x$DEseq2$de
#                                                                                                                                                                                                                                                                              re$diff=ifelse(re$padj<=0.05,ifelse(re$fold>0,"UP","DOWN"),"NO")
#                                                                                                                                                                                                                                                                              return(re)
#                                                                                                                                                                                                                                                                            })
#                                                                                                                                                                                                                                                                            > lapply(res1, function(x)table(x$diff))
#                                                                                                                                                                                                                                                                            [[1]]
#                                                                                                                                                                                                                                                                            < table of extent 0 >
#                                                                                                                                                                                                                                                                              
#                                                                                                                                                                                                                                                                              [[2]]
#                                                                                                                                                                                                                                                                            < table of extent 0 >
#                                                                                                                                                                                                                                                                              
#                                                                                                                                                                                                                                                                              > res1$
#                                                                                                                                                                                                                                                                              res1$
#                                                                                                                                                                                                                                                                              > Dba$contrasts[[1]]$DESeq2
#                                                                                                                                                                                                                                                                            > res1 <- lapply(Dba$contrasts, function(x){
#                                                                                                                                                                                                                                                                              re=x$DESeq2$de
#                                                                                                                                                                                                                                                                              re$diff=ifelse(re$padj<=0.05,ifelse(re$fold>0,"UP","DOWN"),"NO")
#                                                                                                                                                                                                                                                                              return(re)
#                                                                                                                                                                                                                                                                            })
#                                                                                                                                                                                                                                                                            > lapply(res1, function(x)table(x$diff))
#                                                                                                                                                                                                                                                                            [[1]]
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            DOWN     NO     UP 
#                                                                                                                                                                                                                                                                            5926 159038   3645 
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            [[2]]
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            DOWN     NO     UP 
#                                                                                                                                                                                                                                                                            4484 156521   7604 
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            > res2 <- lapply(Dbatest$contrasts, function(x){
#                                                                                                                                                                                                                                                                              re=x$DESeq2$de
#                                                                                                                                                                                                                                                                              re$diff=ifelse(re$padj<=0.05,ifelse(re$fold>0,"UP","DOWN"),"NO")
#                                                                                                                                                                                                                                                                              return(re)
#                                                                                                                                                                                                                                                                            })
#                                                                                                                                                                                                                                                                            > lapply(res1, function(x)table(x$diff))
#                                                                                                                                                                                                                                                                            [[1]]
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            DOWN     NO     UP 
#                                                                                                                                                                                                                                                                            5926 159038   3645 
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            [[2]]
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            DOWN     NO     UP 
#                                                                                                                                                                                                                                                                            4484 156521   7604 
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            > lapply(res2, function(x)table(x$diff))
#                                                                                                                                                                                                                                                                            [[1]]
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            DOWN     NO     UP 
#                                                                                                                                                                                                                                                                            5926 159038   3645 
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            [[2]]
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            DOWN     NO     UP 
#                                                                                                                                                                                                                                                                            4484 156521   7604 
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            [[3]]
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            DOWN     NO     UP 
#                                                                                                                                                                                                                                                                            5402 157199   6008 
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            > tes=dba.report(Dba)
#                                                                                                                                                                                                                                                                            > head(tes
#                                                                                                                                                                                                                                                                                   + )
#                                                                                                                                                                                                                                                                            GRanges object with 6 ranges and 6 metadata columns:
#                                                                                                                                                                                                                                                                              seqnames              ranges strand |      Conc Conc_p30UT Conc_p42UT
#                                                                                                                                                                                                                                                                            <Rle>           <IRanges>  <Rle> | <numeric>  <numeric>  <numeric>
#                                                                                                                                                                                                                                                                              2056     chr1   40300808-40303914      * |   9.63436   10.31464    6.25828
#                                                                                                                                                                                                                                                                            124964     chr5   96940654-96944004      * |   9.16168    6.45043   10.34498
#                                                                                                                                                                                                                                                                            54174    chr14   63797574-63800989      * |   9.10757    9.77289    6.06196
#                                                                                                                                                                                                                                                                            65800    chr16   24017549-24020566      * |   8.13295    5.68357    9.28699
#                                                                                                                                                                                                                                                                            99097     chr2 130770490-130777852      * |  10.00627    8.09692   11.07713
#                                                                                                                                                                                                                                                                            87525    chr19   21437802-21440266      * |   8.61192    6.13071    9.76984
#                                                                                                                                                                                                                                                                            Fold     p-value         FDR
#                                                                                                                                                                                                                                                                            <numeric>   <numeric>   <numeric>
#                                                                                                                                                                                                                                                                              2056   4.01463 1.05680e-68 1.78185e-63
#                                                                                                                                                                                                                                                                            124964  -3.85608 4.53734e-66 3.82516e-61
#                                                                                                                                                                                                                                                                            54174   3.67554 4.99612e-60 2.80795e-55
#                                                                                                                                                                                                                                                                            65800  -3.55563 1.20364e-48 5.07359e-44
#                                                                                                                                                                                                                                                                            99097  -2.93046 2.42299e-48 8.17071e-44
#                                                                                                                                                                                                                                                                            87525  -3.59073 1.61185e-47 4.52951e-43
#                                                                                                                                                                                                                                                                            -------
#                                                                                                                                                                                                                                                                              seqinfo: 24 sequences from an unspecified genome; no seqlengths
#                                                                                                                                                                                                                                                                            > leng
#                                                                                                                                                                                                                                                                            Error: objeto 'leng' no encontrado
#                                                                                                                                                                                                                                                                            > length(tes)
#                                                                                                                                                                                                                                                                            [1] 9571
#                                                                                                                                                                                                                                                                            > Dba$contrasts$
#                                                                                                                                                                                                                                                                              Dba$contrasts$
#                                                                                                                                                                                                                                                                              > Dba$contrasts[[1]]$
#                                                                                                                                                                                                                                                                              Dba$contrasts[[1]]$contrastType  Dba$contrasts[[1]]$group1
#                                                                                                                                                                                                                                                                            Dba$contrasts[[1]]$contrast      Dba$contrasts[[1]]$group2
#                                                                                                                                                                                                                                                                            Dba$contrasts[[1]]$name1         Dba$contrasts[[1]]$DESeq2
#                                                                                                                                                                                                                                                                            Dba$contrasts[[1]]$name2         
#                                                                                                                                                                                                                                                                            > Dba$contrasts[[1]]$name1
#                                                                                                                                                                                                                                                                            [1] "p30UT"
#                                                                                                                                                                                                                                                                            > Dba$contrasts[[1]]$contrastType
#                                                                                                                                                                                                                                                                            [1] "simple"
#                                                                                                                                                                                                                                                                            > Dba$contrasts[[1]]$group1
#                                                                                                                                                                                                                                                                            EV_LPS_R1   EV_UT_R1 p30_LPS_R1 p30_LPS_R2 p30_LPS_R3  p30_UT_R1  p30_UT_R2 
#                                                                                                                                                                                                                                                                            FALSE      FALSE      FALSE      FALSE      FALSE       TRUE       TRUE 
#                                                                                                                                                                                                                                                                            p30_UT_R3 p42_LPS_R1 p42_LPS_R2  p42_UT_R1  p42_UT_R2 
#                                                                                                                                                                                                                                                                            TRUE      FALSE      FALSE      FALSE      FALSE 
#                                                                                                                                                                                                                                                                            > Dbatest$contrasts[[1]]$
#                                                                                                                                                                                                                                                                              Dbatest$contrasts[[1]]$contrastType  Dbatest$contrasts[[1]]$DESeq2
#                                                                                                                                                                                                                                                                            Dbatest$contrasts[[1]]$contrast      
#                                                                                                                                                                                                                                                                            > Dbatest$contrasts[[1]]$contrast
#                                                                                                                                                                                                                                                                            [1]  0  0  0  1  0 -1
#                                                                                                                                                                                                                                                                            > Dbatest$contrasts[[1]]$
#                                                                                                                                                                                                                                                                              Dbatest$contrasts[[1]]$contrastType  Dbatest$contrasts[[1]]$DESeq2
#                                                                                                                                                                                                                                                                            Dbatest$contrasts[[1]]$contrast      
#                                                                                                                                                                                                                                                                            > Dbatest$contrasts[[1]]$
#                                                                                                                                                                                                                                                                              >  dba.show(Dba, bContrast=T)
#                                                                                                                                                                                                                                                                            Factor  Group Samples Group2 Samples2 DB.DESeq2
#                                                                                                                                                                                                                                                                            1 Condition  p30UT       3  p42UT        2      9571
#                                                                                                                                                                                                                                                                            2 Condition p30LPS       3 p42LPS        2     12088
#                                                                                                                                                                                                                                                                            >  dba.show(Dbatest, bContrast=T)
#                                                                                                                                                                                                                                                                            Factor        Group Intercept Condition_EVUT_vs_EVLPS
#                                                                                                                                                                                                                                                                            1 Coefficient 0,0,0,1,0,-1         0                       0
#                                                                                                                                                                                                                                                                            2 Coefficient 0,0,1,0,-1,0         0                       0
#                                                                                                                                                                                                                                                                            3 Coefficient  0,0,0,0,1,0         0                       0
#                                                                                                                                                                                                                                                                            Condition_p30LPS_vs_EVLPS Condition_p30UT_vs_EVLPS Condition_p42LPS_vs_EVLPS
#                                                                                                                                                                                                                                                                            1                         0                        1                         0
#                                                                                                                                                                                                                                                                            2                         1                        0                        -1
#                                                                                                                                                                                                                                                                            3                         0                        0                         1
#                                                                                                                                                                                                                                                                            Condition_p42UT_vs_EVLPS DB.DESeq2
#                                                                                                                                                                                                                                                                            1                       -1      9571
#                                                                                                                                                                                                                                                                            2                        0     12088
#                                                                                                                                                                                                                                                                            3                        0     11410
#                                                                                                                                                                                                                                                                            > Dba_res=dba.report(Dba,contrast = 1:length(Dba$contrasts))
#                                                                                                                                                                                                                                                                            Error: Can only specify one contrast unless requesting a report-based DBA.
#                                                                                                                                                                                                                                                                            > Dba_res=lapply(1:length(Dba$contrasts),function(x)
#                                                                                                                                                                                                                                                                              dba.report(Dba,contrast = x,DataType = 'DBA_DATA_FRAME')
#                                                                                                                                                                                                                                                                              + )
#                                                                                                                                                                                                                                                                            > length(Dba_res)
#                                                                                                                                                                                                                                                                            [1] 2
#                                                                                                                                                                                                                                                                            > head(Dba_res[[1]])
#                                                                                                                                                                                                                                                                            Chr     Start       End      Conc Conc_p30UT Conc_p42UT      Fold
#                                                                                                                                                                                                                                                                            2056    chr1  40300808  40303914  9.634357  10.314641   6.258283  4.014626
#                                                                                                                                                                                                                                                                            124964  chr5  96940654  96944004  9.161677   6.450431  10.344975 -3.856078
#                                                                                                                                                                                                                                                                            54174  chr14  63797574  63800989  9.107569   9.772895   6.061959  3.675540
#                                                                                                                                                                                                                                                                            65800  chr16  24017549  24020566  8.132948   5.683568   9.286990 -3.555629
#                                                                                                                                                                                                                                                                            99097   chr2 130770490 130777852 10.006270   8.096923  11.077127 -2.930462
#                                                                                                                                                                                                                                                                            87525  chr19  21437802  21440266  8.611920   6.130711   9.769843 -3.590726
#                                                                                                                                                                                                                                                                            p-value          FDR
#                                                                                                                                                                                                                                                                            2056   1.056798e-68 1.781846e-63
#                                                                                                                                                                                                                                                                            124964 4.537337e-66 3.825157e-61
#                                                                                                                                                                                                                                                                            54174  4.996123e-60 2.807955e-55
#                                                                                                                                                                                                                                                                            65800  1.203643e-48 5.073595e-44
#                                                                                                                                                                                                                                                                            99097  2.422989e-48 8.170706e-44
#                                                                                                                                                                                                                                                                            87525  1.611848e-47 4.529509e-43
#                                                                                                                                                                                                                                                                            > saveRDS(Dbatest,"Dba_first_3comps.RDS")
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            > 
#                                                                                                                                                                                                                                                                              > remaining_ctsts=list(c(0,0,1,0,0,0),
#                                                                                                                                                                                                                                                                                                     c(0,-1,0,0,0,1),
#                                                                                                                                                                                                                                                                                                     c(0,-1,0,1,0,0),
#                                                                                                                                                                                                                                                                                                     c(0,0,1,-1,0,0),
#                                                                                                                                                                                                                                                                                                     c(0,0,0,0,1,-1),
#                                                                                                                                                                                                                                                                                                     c(0,-1,0,0,0,0))
#                                                                                                                                                                                                                                                                            > for (ct in remaining_ctsts) {
#                                                                                                                                                                                                                                                                              newdba=dba.contrast(newdba,contrast = ct)
#                                                                                                                                                                                                                                                                            }
#                                                                                                                                                                                                                                                                            Error in dba.contrast(newdba, contrast = ct) : 
#                                                                                                                                                                                                                                                                              objeto 'newdba' no encontrado
#                                                                                                                                                                                                                                                                            > newdba=DBA
#                                                                                                                                                                                                                                                                            > DBA$DESeq2
#                                                                                                                                                                                                                                                                            $names
#                                                                                                                                                                                                                                                                            [1] "Intercept"                 "Condition_EVUT_vs_EVLPS"  
#                                                                                                                                                                                                                                                                            [3] "Condition_p30LPS_vs_EVLPS" "Condition_p30UT_vs_EVLPS" 
#                                                                                                                                                                                                                                                                            [5] "Condition_p42LPS_vs_EVLPS" "Condition_p42UT_vs_EVLPS" 
#                                                                                                                                                                                                                                                                            
#                                                                                                                                                                                                                                                                            > for (ct in remaining_ctsts) {
#                                                                                                                                                                                                                                                                              newdba=dba.contrast(newdba,contrast = ct)
#                                                                                                                                                                                                                                                                            }
#                                                                                                                                                                                                                                                                            > length(newdba$contrasts)
#                                                                                                                                                                                                                                                                            [1] 6
#                                                                                                                                                                                                                                                                            > newdba=dba.analyze(newdba)
#                                                                                                                                                                                                                                                                            Applying Blacklist/Greylists...
#                                                                                                                                                                                                                                                                            No genome detected.
#                                                                                                                                                                                                                                                                            Normalize DESeq2 with defaults...
#                                                                                                                                                                                                                                                                            Analyzing...
#                                                                                                                                                                                                                                                                            gene-wise dispersion estimates
#                                                                                                                                                                                                                                                                            mean-dispersion relationship
#                                                                                                                                                                                                                                                                            final dispersion estimates
#                                                                                                                                                                                                                                                                            > getw()
#                                                                                                                                                                                                                                                                            Error in getw() : no se pudo encontrar la función "getw"
#                                                                                                                                                                                                                                                                            > getwd()
#                                                                                                                                                                                                                                                                            [1] "/mnt/beegfs/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind"
#                                                                                                                                                                                                                                                                            > list.files("/tmp/llorenzi/1962276")
#                                                                                                                                                                                                                                                                            [1] "RtmpFyVIaw"
#                                                                                                                                                                                                                                                                            > list.files()
#                                                                                                                                                                                                                                                                            [1] "CEBPa_DiffBind.nosummarizeoverlap.csuc.keptlength.R"
#                                                                                                                                                                                                                                                                            [2] "Dba_first_3comps.RDS"                               
#                                                                                                                                                                                                                                                                            [3] "DBA.keptlength.RDS"                                 
#                                                                                                                                                                                                                                                                            [4] "PCA_CEBPa_minOverlap_1.keptlength.pdf"              
#                                                                                                                                                                                                                                                                            [5] "plotMdba_minOverlap_1.pdf"                          
#                                                                                                                                                                                                                                                                            [6] "sample_sheet_csuc.csv"                              
#                                                                                                                                                                                                                                                                            > saveRDS(newdba,"Dba_remainingcomps.RDS")
#                                                                                                                                                                                                                                                                            > 
#                                                                                                                                                                                                                                                                              
#                                                                                                                                                                                                                                                                              