
#bws is a vector with bigwig paths, peaks is a Granges object
bws=c("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/ER/bigWig_files/nov21/p30_ER.sorted.bam.CPM.bw",
      "~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/ER/bigWig_files/nov21/p42_ER.sorted.bam.CPM.bw")
p30_peaks=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/ER/p30_vs_p42/bed_files/p30_ER_UT_peaks.bed")


peaks=p30_peaks

#introduce false peaks just to test what happens if missing coordinates
peaks=peaks[1:100,]
peaks=rbind(peaks,c("chr2",181755117,181756517,"chr2:181755117-181756517",8.234545,"."))
#this peak is entirely out of coordinates

chr1_size=195154279


#add a peak that is partially out of coordinates
peaks=rbind(peaks,c("chr1",chr1_size - 1000,chr1_size +1000,
                    paste0("chr1:",chr1_size - 1000,"-",
                           chr1_size+1000),4.345323,"."))


#add this peak that I saw that failed before
ptk="chr1_GL456239v1_random:38335-38511"
peaks=rbind(peaks,p30_peaks[p30_peaks$V4==ptk,])

region_to_plot=c(-1000,1000)
##
#lower limit
llimit=region_to_plot[1]
#upper limit
ulimit=region_to_plot[2]
#region width
rwidth=length(llimit:ulimit)
chrs=peaks$V1
midpoints=midpoint(peaks$V2,peaks$V3)
llimit=1000
rwidth

#make peaks that have wrong coordinates to test import function
peaks
peaks_gr=granges_from_midpoints(chrs,midpoints,llimit,rwidth)

bw=bws[2]
rm_overlap=F
rm_incomplete=T

rm_overlap=T
rm_incomplete=T

rm_overlap=T
rm_incomplete=F


library(rtracklayer)  
rwidth=width(extgr)[1]
extgr=peaks_gr

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