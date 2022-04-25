## ---------------------------
##
##
## Purpose of script: process DiffBind results from hpc
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-03-22
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes:#I discovered partially the origin discrepancies between differential results 
#as computed by DiffBind dba.show() and those obtained after retrieving
#the DESeq2 object from DiffBind () and then running results() 
#(see these posts: https://support.bioconductor.org/p/9135706/ ; 
#https://support.bioconductor.org/p/9135779/#9135954)
#The thing is that when you run DESeq2 results() function
#with different alpha values, the FDR changes slightly for all peaks
## If you get the results for all peaks with diffbind using dba.report with th=1,
## the results are almost equivalent to those retrieved with results with 0.1 (default)
##
## ---------------------------

options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
require(DiffBind)
## ---------------------------
setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/DiffBind/DESeq2/")
Dba3first=readRDS("Dba_first_3comps.RDS")
Dbaremaining <- readRDS("Dba_remainingcomps.RDS")

#retrieve peakset info (coordinates and counts)
peakset=dba.peakset(Dba3first,bRetrieve = T,
                    DataType = DBA_DATA_FRAME)
peakset$ID <- paste0(peakset$CHR,":",peakset$START,"-",peakset$END)
#retrieve DESeq2 dataset
dds <- dba.analyze(Dba3first, bRetrieveAnalysis=DBA_DESEQ2)

#generate differential results tables for all contrasts
coeffs=Dba3first$DESeq2$names
res_df <- list()
diffbind_res_df <- list()
ctsts <- list()
for(x in 1:length(Dba3first$contrasts)){
  ct=Dba3first$contrasts[[x]]$contrast
  print(ct)
  comp=coeffs[ct!=0][order(ct[ct!=0],decreasing = T)]
  cw=sort(ct[ct!=0],decreasing = T)
  comp[cw<0] <- paste0("-",comp[cw<0])
  
  comp <- paste(comp,collapse = "")
  ctsts[[comp]] <- ct
  diffbind_res_df[[comp]] <- as.data.frame(dba.report(Dba3first,contrast = x,th = 1))
  res_df[[comp]] <- as.data.frame(results(dds2,contrast = ct))
}

for(x in 1:length(Dbaremaining$contrasts)){
  ct=Dbaremaining$contrasts[[x]]$contrast
  print(ct)
  comp=coeffs[ct!=0][order(ct[ct!=0],decreasing = T)]
  cw=sort(ct[ct!=0],decreasing = T)
  comp[cw<0] <- paste0("-",comp[cw<0])
  
  comp <- paste(comp,collapse = "")
  ctsts[[comp]] <- ct
  diffbind_res_df[[comp]] <- as.data.frame(dba.report(Dbaremaining,contrast = x,th = 1))
  res_df[[comp]] <- as.data.frame(results(dds2,contrast = ct))
}

table(res_df$`Condition_p30UT_vs_EVLPS-Condition_p42UT_vs_EVLPS`$padj<=0.05)
table(diffbind_res_df$`Condition_p30UT_vs_EVLPS-Condition_p42UT_vs_EVLPS`$FDR<=0.05)

res_df <- lapply(res_df, function(x){rownames(x)=peakset$ID
return(x)})
View(res_df$`Condition_p30UT_vs_EVLPS-Condition_p42UT_vs_EVLPS`)

diffbind_res_df <- lapply(diffbind_res_df,function(x){
  x$ID=paste0(x$seqnames,":",x$start,"-",x$end)
  colnames(x)[6:9] <- paste0(colnames(x)[6:9],".DiffBind")
  return(x)
})

res_all <- lapply(names(diffbind_res_df),function(na) {
  deseq=res_df[[na]]
  dbind=diffbind_res_df[[na]]
  dbind <- dbind[match(rownames(deseq),dbind$ID),]
  tmp=cbind(dbind[,c(10,1:5)],deseq,dbind[6:9])
  return(tmp)
  
})
cnames=names(diffbind_res_df)
cnames[str_count(cnames,pattern = "_vs_")==2] <- gsub("-","_vs_",gsub("_vs_EVLPS","",cnames[str_count(cnames,pattern = "_vs_")==2]))
cnames <- gsub("Condition_","",cnames)
cnames[grep("-",cnames)] <- sapply(cnames[grep("-",cnames)],
                                   function(x){
                                     paste(gsub("-","",
                                                rev(strsplit(x,split = "_vs_")[[1]])),
                                           collapse = "_vs_")})
names(res_all) <- cnames

View(res_all$p30UT_vs_p42UT)

res_all <- lapply(res_all,function(x){
  x$diffDESeq2=ifelse(x$padj<=0.05&!is.na(x$padj),ifelse(x$log2FoldChange>0,"UP","DOWN"),"NO")
  x$diffDiffBind=ifelse(x$FDR<=0.05&!is.na(x$FDR),ifelse(x$Fold.DiffBind>0,"UP","DOWN"),"NO")
  return(x)
})

for (na in names(res_all)) {
  nam=gsub("_","",na)
  dir.create(nam)
  res=res_all[[na]]
  capture.output(print("diffDESeq2 vs diffDiffBind"),
                 file=paste0(nam,"/",nam,"_summary_res.txt"))
  
  capture.output(table(res$diffDESeq2,res$diffDiffBind),
                 file=paste0(nam,"/",nam,"_summary_res.txt"),append = T)
  res <- res[order(res$padj),]
  write.csv(res,paste0(nam,"/",nam,"_DESeq2_DiffBind_res.csv"))
}
