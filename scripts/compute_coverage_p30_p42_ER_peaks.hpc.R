## ---------------------------
##
##
## Purpose of script: compute coverage objects for p30 and p42 ER peaks
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-05-16
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------
cargs=commandArgs(trailingOnly = T)
bw_paths=cargs[1] # /scratch/llorenzi/coverage_p30_p42/list_bigWig_files.txt  
peaksobj=cargs[2] # /scratch/llorenzi/coverage_p30_p42/peak_sets_GR.RDS 

options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
## ---------------------------
require(rtracklayer)

#functions
midpoint <- function(starts,ends){
  mp=sapply(1:length(starts),function(i){
    s=as.numeric(starts[i])
    e=as.numeric(ends[i])
    wid=e - s +1
    return(s + round(wid/2)-1)
  })
  return(mp)
}

granges_from_midpoints <-  function(chrs,midpoints,llimit,rwidth){
  library(GenomicRanges)
  return(GRanges(chrs,
                 IRanges(start = midpoints + llimit,
                         width = rwidth)))
} 

compute_coverage_around_peaks=function(bws,extgr,rm_overlap=F,rm_incomplete=T){
  library(rtracklayer)  
  rwidth=width(extgr)[1]
  #things to consider: overlap between peaks report on this
  #number of peaks: if too long then split
  #are peaks complete after importing?
  
  #calculate and report on overlaps
  print(paste0("There are ",length(extgr)," sequences to process"))
  print(paste0("Regions width is ",rwidth))
  ###Remove possible overlapping regions 
  tfo=findOverlaps(extgr,extgr)
  overlaps=data.frame(q=queryHits(tfo),s=subjectHits(tfo))
  overlaps <- overlaps[overlaps$q!=overlaps$s,]
  #first remove the repeated rows, keep only those in which q is smaller than s 
  overlapsf <- overlaps[overlaps$q<overlaps$s,]
  rows_to_remove <- unique(overlapsf$s)
  if(length(rows_to_remove)>0){
    print(paste0("There are ",length(rows_to_remove)," overlaps"))
    if(rm_overlap){extgr <- extgr[-rows_to_remove]
    print("For overlaps, one of the peaks will be removed")
    }
  }
  
  # 
  chunk_length=round(10005000/rwidth)
  chunks=list()
  if(ceiling(length(extgr)/chunk_length)>1){
    chunks=ceiling(length(extgr)/chunk_length)
    cl=ceiling(length(extgr)/chunks)
    
    i=1
    j=cl
    chunks=list()
    n=1
    while (j<length(extgr)) {
      chunks[[n]]=i:j
      i=j+1
      j=i+cl -1
      n=n+1
    }
    chunks[[n]]=i:length(extgr)
  }else chunks[[1]]=1:length(extgr)
  
  print(paste0("Regions will be processed in ",length(chunks)," chunks"))
  
  all_TSSposCPM=list()
  for (cn in 1:length(chunks)) {
    print(paste0("Processing chunk ",cn," of ", length(chunks)))
    ext_gr=extgr[chunks[[cn]]]
    
    #1) generate id per position in origin
    posIDorigin <- paste0(as.character(rep(seqnames(ext_gr),each=rwidth)),"_",
                          unlist(apply(as.data.frame(ranges(ext_gr)),1,function(x)return(x[1]:x[2]))))
    #2)generate corresponding peakIDs for each position in the original peaks:      
    chrs_origin <- as.character(rep(seqnames(ext_gr),each=rwidth))
    ranges_origin <- rep(unlist(apply(as.data.frame(ranges(ext_gr)),
                                      1,function(x)return(paste(x[1:2],collapse = "-")))),
                         each=rwidth)
    peakIDorigin <- paste0(chrs_origin,":",ranges_origin)
    
    if(!rm_overlap){
      #3) check peaks with overlap
      print("Duplicated positions because of regions overlap")
      print(table(duplicated(posIDorigin)))
      #peaks with duplications
      dups=posIDorigin[duplicated(posIDorigin)]
      print(paste0("there are ", length(unique(dups))," duplicated positions"))
      print("times duplicated:")
      print(table(table(dups)))
      peakIDdups=peakIDorigin[posIDorigin%in%dups]
      posIDdups=posIDorigin[posIDorigin%in%dups]
      td=table(peakIDorigin[posIDorigin%in%dups])
      tdt=table(posIDdups,peakIDdups)
      tdt=tdt[!duplicated(tdt)]  
      print(paste0("There are ",length(td)," peaks with overlap"))
      
    }
    
    
    for (bw in bws) {
      print(paste0("Importing coverage from ", bw, " for chunk ",cn))
      tmp_tssCPM=import(bw,which=ext_gr)
      #extend the GRanges object so each position around -Xkb - +X kb of each TSS has a CPM value
      tmp_all_TSSposCPM <- rep(tmp_tssCPM$score,width(ranges(tmp_tssCPM)))
      print(paste0("Finished import of ", bw, " for chunk ",cn))
      print("Checking if incomplete peaks")
      #####Check that each feature have all positions. If incomplete, remove!####
      
      #generate single-position ids in the coverage object
      posIDCovobj <- paste0(as.character(rep(seqnames(tmp_tssCPM),width(ranges(tmp_tssCPM)))),"_",
                            unlist(apply(as.data.frame(ranges(tmp_tssCPM)),1,function(x)return(x[1]:x[2]))))
      
      #4) Use posIDorigin to select the cpm score of each position of each peak of origin 
      #(this conveniently reorders the cpm vector to match original positions):
      #If for some reason there are missing positions (unmatched origin positions) 
      #these will be removed in next step
      tmp_all_TSSposCPM <- tmp_all_TSSposCPM[match(posIDorigin,posIDCovobj)]
      
      #are there any missing positions?
      missing_positions=sum(is.na(tmp_all_TSSposCPM))
      if(missing_positions>0){
        print(paste0("There are ",missing_positions," missing positions"))
        print(paste0("Missing positions:"))
        print(posIDorigin[is.na(tmp_all_TSSposCPM)])
        ti=table(peakIDorigin[is.na(tmp_all_TSSposCPM)])
        print(paste0(length(ti)," incomplete region(s) found"))
        print("number of incomplete positions per region: ")
        print(ti)
        print("This is most likely because these regions are not present in the bigWig file")
        
        if(rm_incomplete){
          print("Incomplete regions will be removed")
          peaks_to_remove=names(ti)
          tmp_all_TSSposCPM=tmp_all_TSSposCPM[!peakIDorigin%in%peaks_to_remove]
          
        }else {
          print("Zeros will be added to these positions")
          tmp_all_TSSposCPM[is.na(tmp_all_TSSposCPM)] <- 0
        }
        
      }   else print("There are no missing positions")
      
      print("There are")
      print(length(tmp_all_TSSposCPM)/rwidth)
      print("regions to plot")
      
      
      all_TSSposCPM[[basename(bw)]]=c(all_TSSposCPM[[basename(bw)]],tmp_all_TSSposCPM)
      
    }
    print(paste0("Finished data import chunk ",cn, " of ", length(chunks)))
    
  }
  return(all_TSSposCPM)
}

plot_densities <- function(dfin,lwd=2,main="",
                           xlab="distance from Peak center",
                           ylab="mean CPM",xvals=llimit:ulimit,
                           cols=cbPalette[6:7],
                           sampleTable=sample_table,
                           groupCol=group_feat){
  #Define sample sets: has to be a factor of two groups:
  #group_feat="LPS" #to test
  sampleTable[,groupCol] <- as.factor(sampleTable[,groupCol])
  sgr <- levels(sampleTable[,groupCol])
  vcols <- cols[as.numeric(as.factor(sampleTable[,groupCol]))]
  
  #we use it to generate average vectors for each position*peak
  #across all samples from each group
  
  maxyval=max(apply(dfin, 2,max))
  minyval=min(apply(dfin, 2,min))
  plot(xvals, dfin[,1],col=vcols[1],
       type="s",ylim=c(minyval,maxyval+0.02),ylab=ylab,
       xlab=xlab,lwd=lwd,main=main
  )
  
  legend("topright",lwd = lwd,legend = sgr,
         title = groupCol,bty = "n",col = cols)
  for (sa in 2:ncol(dfin)) {
    points(xvals, dfin[,sa],
           type="s",col=vcols[sa],lwd=2)
  }
  
}

region_to_plot=c(-3000,3000)
##
#lower limit
llimit=region_to_plot[1]
#upper limit
ulimit=region_to_plot[2]
#region width
rwidth=length(llimit:ulimit)

bws=read.table(bw_paths)[,1]
list_peak_sets_gr=readRDS(peaksobj)

for (bw in bws) {
  bw_tmp_CPMs=lapply(list_peak_sets_gr,function(x)compute_coverage_around_peaks(bws = bw,extgr = x,rm_overlap = F,rm_incomplete = T))
  saveRDS(bw_tmp_CPMs,paste0(basename(bw),".covObj.RDS"))
}