tmp_tssCPM_numlist=import(bw,which=extended_granges,as="NumericList")
tmp_tssCPM
length(extended_granges)
lapply(tmp_tssCPM,length)
unlist(tmp_tssCPM_numlist) == all_TSSposCPM
  
names(extended_granges)
#1) generate single-position ids in the peaks object:
posIDorigin <- paste0(as.character(rep(seqnames(extended_granges),each=rwidth)),"_",
                      unlist(apply(as.data.frame(ranges(extended_granges)),1,function(x)return(x[1]:x[2]))))
#2)generate corresponding peakIDs for each position in the original peaks:
chrs_origin <- as.character(rep(seqnames(extended_granges),each=rwidth))
ranges_origin <- rep(unlist(apply(as.data.frame(ranges(extended_granges)),
                                  1,function(x)return(paste(x[1:2],collapse = "-")))),
                     each=rwidth)
peakIDorigin <- paste0(chrs_origin,":",ranges_origin)


#3) In the same way, generate single-position ids in the coverage object
posIDCovobj <- paste0(as.character(rep(seqnames(tmp_tssCPM),width(ranges(tmp_tssCPM)))),"_",
                      unlist(apply(as.data.frame(ranges(tmp_tssCPM)),1,function(x)return(x[1]:x[2]))))
#4) Use both to assign peak of origin for each position in the coverage object
peakIDCovobj <- peakIDorigin[match(posIDCovobj,posIDorigin)]


#assign names to tmp_tssCPM_numlist
numlistmetadata=as.data.frame(metadata(tmp_tssCPM_numlist))
names(tmp_tssCPM_numlist)=paste0(numlistmetadata$ranges.seqnames,":",
                                 numlistmetadata$ranges.start,"-",
                                 numlistmetadata$ranges.end)

table(names(tmp_tssCPM_numlist)%in%peakIDorigin)
table(names(tmp_tssCPM_numlist)==peakIDorigin)
peakIDnumlist=rep(names(tmp_tssCPM_numlist),each=rwidth)
posIDnumlist=paste0(as.character(rep(numlistmetadata$ranges.seqnames,numlistmetadata$ranges.width)),"_",
                    unlist(apply(numlistmetadata[,2:3],1,function(x)return(x[1]:x[2]))))

table(posIDnumlist%in%posIDorigin)

unlistednumlist=unlist(tmp_tssCPM_numlist)
tmp=names(unlistednumlist)
table(tmp==peakIDnumlist)
names(unlistednumlist)=posIDnumlist
all_TSSposCPM_reord=all_TSSposCPM[match(posIDnumlist,names(all_TSSposCPM))]
table(unlistednumlist==all_TSSposCPM[match(names(unlistednumlist),
                                           names(all_TSSposCPM))])
table(all_TSSposCPM_reord==unlistednumlist)
all_TSSposCPM_reord[all_TSSposCPM_reord!=unlistednumlist]
tmp_tssCPM_numlist[["chr1:133621070-133623070"]]
