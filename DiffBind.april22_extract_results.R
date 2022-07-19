## ---------------------------
##
##
## Purpose of script: process diffbind resutls
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-04-26
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
require(DiffBind)
## ---------------------------
setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/DiffBind/new_apr22")

Dba=readRDS("Mdba.analysis.apr22.RDS")
deseq2_resl <- list()
edger_resl <- list()
ctsts=list("p30UT_vs_p42UT"=c(0,0,0,1,0,-1),   
           "p30LPS_vs_p42LPS"=c(0,0,1,0,-1,0),
           "p42LPS_vs_EVLPS"=c(0,0,0,0,1,0),   
           "p30LPS_vs_EVLPS"=c(0,0,1,0,0,0),  
           "p42UT_vs_EVUT"=c(0,-1,0,0,0,1),    
           "p30UT_vs_EVUT"=c(0,-1,0,1,0,0),   
           "p30LPS_vs_p30UT"=c(0,0,1,-1,0,0),  
           "p42LPS_vs_p42UT"=c(0,0,0,0,1,-1), 
           "EVLPS_vs_EVUT"=c(0,-1,0,0,0,0) )

table(unlist(lapply(Dba$contrasts,function(x)x$contrast))==unlist(ctsts))

for(x in 1:length(Dba$contrasts)){
  ct=Dba$contrasts[[x]]$contrast
  print(ct)
  comp= names(ctsts)[x]
  deseq2_resl[[comp]] <- as.data.frame(dba.report(Dba,contrast = x,method=DBA_DESEQ2, th = 1))
  edger_resl[[comp]] <- as.data.frame(dba.report(Dba,contrast = x,method=DBA_EDGER, th = 1))
}

ress=list(DESeq2=deseq2_resl,edgeR=edger_resl)
ress <- lapply(ress,function(y)lapply(y,function(x){
  x$ID=paste0(x$seqnames,":",x$start,"-",x$end)
  x$diff=ifelse(x$FDR<=0.05&!is.na(x$FDR),ifelse(x$Fold>0,"UP","DOWN"),"NO")
  return(x)
}))



for (nam in names(ress[[1]])) {
  
  dir.create(nam)
  for (met in names(ress)) {
    capture.output(print(paste0(met," differential peaks")),
                   file=paste0(nam,"/",nam,"_",met,"_summary_res.txt"))
    
    capture.output(table(ress[[met]][[nam]]$diff),
                   file=paste0(nam,"/",nam,"_",met,"_summary_res.txt"),append = T)
    
    res = ress[[met]][[nam]]
    res <- res[order(res$FDR),]
    write.csv(res,paste0(nam,"/",nam,"_",met,"_DiffBind_res.csv"))
    
  }
}  
  