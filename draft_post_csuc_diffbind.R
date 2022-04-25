Dba3first=readRDS("Dba_first_3comps.RDS")

coeffs=Dba3first$DESeq2$names
res_df <- list()
ctsts <- list()
for(x in 1:length(Dba3first$contrasts)){
  ct=Dba3first$contrasts[[x]]$contrast
  print(ct)
  comp=coeffs[ct!=0][order(ct[ct!=0],decreasing = T)]
  cw=sort(ct[ct!=0],decreasing = T)
  comp[cw<0] <- paste0("-",comp[cw<0])

  comp <- paste(comp,collapse = "")
  ctsts[[comp]] <- ct
  res_df[[comp]] <- as.data.frame(dba.report(Dba3first,contrast = x,th = 1))
}

View(res_df$`Condition_p30UT_vs_EVLPS-Condition_p42UT_vs_EVLPS`)

Dbaremaining <- readRDS("Dba_remainingcomps.RDS")
for(x in 1:length(Dbaremaining$contrasts)){
  ct=Dbaremaining$contrasts[[x]]$contrast
  print(ct)
  comp=coeffs[ct!=0][order(ct[ct!=0],decreasing = T)]
  cw=sort(ct[ct!=0],decreasing = T)
  comp[cw<0] <- paste0("-",comp[cw<0])
  
  comp <- paste(comp,collapse = "")
  ctsts[[comp]] <- ct
  res_df[[comp]] <- as.data.frame(dba.report(Dbaremaining,contrast = x,th = 1))
}

cnames=names(res_df)
cnames <- gsub("Condition_","",cnames)
cnames[grep("-",cnames)] <- gsub("-","_vs_",gsub("_vs_EVLPS","",cnames[grep("-",cnames)]))

ndf=lapply(res_df, function(x)x$diff=ifelse(x$FDR<=0.05,ifelse(x$fo)))

te$diff <- ifelse(te$FDR<=0.05)
cnames[str_count(cnames,pattern = "_vs_")==2] <- gsub("-","_vs_",gsub("_vs_EVLPS","",cnames[str_count(cnames,pattern = "_vs_")==2]))
cnames[grep("-",cnames)] <- sapply(cnames[grep("-",cnames)],
                                   function(x){
                                     paste(gsub("-","",
                                                rev(strsplit(x,split = "_vs_")[[1]])),
                                           collapse = "_vs_")})

ctsts
dds <- dba.analyze(Dba3first, bRetrieveAnalysis=DBA_DESEQ2)
library(DESeq2)
results(dds,alpha = 0.05)
summary(results(dds, alpha = 0.05))
te=res_df$`-Condition_EVUT_vs_EVLPS`

te$diff <- ifelse(te$FDR<0.05,ifelse(te$Fold>0,"UP","DOWN"),"NO")
table(te$diff)

dds2 <- dba.analyze(Dbaremaining,bRetrieveAnalysis = DBA_DESEQ2)

summary(results(dds2, alpha = 0.05,contrast = ctsts[[1]]))
summary(results(dds2, alpha = 0.05,contrast = ctsts[[2]]))
summary(results(dds2, independentFiltering=F,contrast = ctsts[[9]]))
summary(results(dds2, independentFiltering=F,contrast = c("Condition","EVLPS","EVUT")))
summary(results(dds2, contrast = c("Condition","EVLPS","EVUT")))
res=as.data.frame(results(dds2, alpha = 0.05,contrast = ctsts[[1]]))
anyNA(res)
dba.show(Dba3first, bContrasts=TRUE)
dba.show(Dba3first, bContrasts=TRUE,th=0.1)
dba.show(Dbaremaining, bContrasts=TRUE)
dba.show(Dbaremaining,bContrasts = T,th=0.1)
summary(results(dds2, independentFiltering = T,contrast = c("Condition","EVLPS","EVUT")))
dba.show(Dbaremaining,bContrasts = T,th=0.05)
summary(results(dds2, independentFiltering = F,contrast = c("Condition","EVLPS","EVUT")))

sapply(ctsts, function(x)summary(results(dds2,contrast = x)))
peakset=dba.peakset(Dba3first,bRetrieve = T,DataType = DBA_DATA_FRAME)

peakset$ID=paste0(peakset$CHR,":",peakset$START,"-",peakset$END)
rownames(res)=peakset$ID
resdb=res_df$`Condition_p30UT_vs_EVLPS-Condition_p42UT_vs_EVLPS`
resdb$ID=paste0(resdb$seqnames,":",resdb$start,"-",resdb$end)
table(resdb$ID%in%rownames(res))
View(res[is.na(res$padj),])
naids <- rownames(res[is.na(res$padj),])
View(resdb[match(naids,resdb$ID),])

table(resdb$FDR<=0.05)
table(res$padj<=0.05)
table(resdb$FDR<=0.1)
table(res$padj<=0.1)


summary(results(dds2, independentFiltering = T,contrast = ctsts[[1]]))
summary(results(dds2, contrast = ctsts[[1]]))

summary(results(dds2, alpha = 0.05, contrast = ctsts[[1]]))
res2=as.data.frame(results(dds2,contrast = ctsts[[1]]))
rownames(res2)=peakset$ID
table(res$padj<=0.1)

table(res2$padj<=0.1)
res3=as.data.frame(results(dds2,alpha = 0.99,contrast = ctsts[[1]]))
table(res3$padj<=0.05)
table(resdb$FDR<=0.05)
table(resdb$FDR<=0.1)
table(res3$padj<=0.1)
rownames(res3) <- rownames(res)
View(resdb[match(rownames(res3),resdb$ID),])

dbids=resdb$ID[resdb$FDR<=0.05]
dsids=rownames(res[res$padj<=0.05&!is.na(res$padj),])
table(dbids%in%dsids)
table(dsids%in%dbids)

View(resdb[match(dsids[!dsids%in%dbids],resdb$ID),])
View(res[rownames(res)%in%dsids[!dsids%in%dbids],])
View(res3[rownames(res)%in%dsids[!dsids%in%dbids],])
View(res2[rownames(res)%in%dsids[!dsids%in%dbids],])
resdb0.05=as.data.frame(dba.report(Dba3first,contrast = 1))
resdb0.05$ID=paste0(resdb0.05$seqnames,":",resdb0.05$start,"-",
                    resdb0.05$end)
View(resdb0.05[match(dsids[dsids%in%dbids],resdb0.05$ID),])
View(resdb[match(dsids[dsids%in%dbids],resdb$ID),])

View(res[rownames(res)%in%dsids[dsids%in%dbids],])

#I discovered partially the origin discrepancies between differential results 
#as computed by DiffBind dba.show() and those obtained after retrieving
#the DESeq2 object from DiffBind () and then running results() 
#(see these posts: https://support.bioconductor.org/p/9135706/ ; 
#https://support.bioconductor.org/p/9135779/#9135954)
#The thing is that when you run DESeq2 results() function
#with different alpha values, the FDR changes slightly for all peaks
#so, for example lets take a peak that I found that changes:
res=as.data.frame(results(dds2, alpha = 0.05,contrast = ctsts[[1]]))
res2=as.data.frame(results(dds2,contrast = ctsts[[1]])) #default alpha value is 0.1
res3=as.data.frame(results(dds2,alpha = 0.99,contrast = ctsts[[1]]))
table(res$padj<=0.05)
# FALSE   TRUE 
# 149073   9728 
table(res2$padj<=0.05)
# FALSE   TRUE 
# 159037   9571 
table(res3$padj<=0.05)
# FALSE   TRUE 
# 159037   9571 
table(res2$padj==res3$padj)
table(res$padj==res2$padj)
res0.5=as.data.frame(results(dds2, alpha = 0.01,contrast = ctsts[[1]]))
table(res0.5$padj<=0.05)
table(res$padj==res0.5$padj)
res1.5=as.data.frame(results(dds2, alpha = 0.075,contrast = ctsts[[1]]))
table(res1.5$padj<=0.05)
table(res1.5$padj==res2$padj)
res0.055=as.data.frame(results(dds2, alpha = 0.055,contrast = ctsts[[1]]))
table(res0.055$padj<=0.05)
table(res0.055$padj==res2$padj)
table(res0.055$padj==res$padj)

pa=res0.055$padj
co=0.055

resall=as.data.frame(results(dds2, alpha = 0.01,contrast = ctsts[[1]]))

while (all(pa==res$padj,na.rm = T)) {
  co=co+0.005 
  restmp=as.data.frame(results(dds2, alpha = co,contrast = ctsts[[1]]))
  pa=restmp$padj
}
res0.06=as.data.frame(results(dds2, alpha = 0.06,contrast = ctsts[[1]]))
table(res0.06$padj<=0.05)

