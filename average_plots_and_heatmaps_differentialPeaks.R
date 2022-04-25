## ---------------------------
##
##
## Purpose of script:average y heatmaps de differential ChIP-seq peaks CEBPa
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-03-17
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
library(rtracklayer)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## ---------------------------
peaksetsdir="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/DESeq2_results/"
metadata=read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/H3K27ac_sample_metadata.txt")
resL <- list.files(peaksetsdir,pattern = "DESeq2_res.noNAs.csv")

bigwigfilesdir="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/bigWig_files/"
list.files(bigwigfilesdir)
bigWigFiles <- list.files(bigwigfilesdir,full.names = T)
bigWigFiles <- bigWigFiles[grep("INP",bigWigFiles,invert = T)]
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/DESeq2_results/ChIP-seq_coverage_diffPeaks_plots/")
padjcutoff=0.05

region_to_plot=c(-1000,1000)
##
#lower limit
llimit=region_to_plot[1]
#upper limit
ulimit=region_to_plot[2]
#region width
rwidth=length(llimit:ulimit)
for (comp in resL) {
  print(comp)
  #parse metadata
  grps <-  unlist(strsplit(gsub("_DESeq2_res.noNAs.csv","",comp),split = "vs"))
  samples <- lapply(grps, function(x)return(metadata$samples[metadata$group==x]))
  names(samples) <- grps
  bwfiles=unlist(lapply(samples, function(x)sapply(x, function(y)grep(y,bigWigFiles,value = T))))
  grp <- rep(names(samples),lapply(samples, length))
  
  #generate peaks midpoints and GRanges
  deseq2res=read.csv(paste0(peaksetsdir,"/",comp))
  upPeaks <- deseq2res %>% filter(padj<=padjcutoff , !is.na(padj) ,log2FoldChange>0)
  downPeaks <- deseq2res %>% filter(padj<=padjcutoff , !is.na(padj) ,log2FoldChange<0)
  
  pl <- list(UPpeaks=upPeaks,
             DOWNpeaks=downPeaks)
  pl <- pl[sapply(pl, function(x)nrow(x)!=0)]
  erl <- sapply(names(pl),function(nam) {
    pp=pl[[nam]]
    chrs <- sapply(strsplit(pp$X,":"),function(x)x[1])
    coords <- strsplit(sapply(strsplit(pp$X,":"),function(x)x[2]),
                       split = "-")
    midpoints <- sapply(coords, function(x){
      x=as.numeric(x)
      wid=x[2] - x[1]+1
      return(x[1] + round(wid/2)-1)
    }
    )
    return( GRanges(chrs,
                                IRanges(start = midpoints + llimit,
                                        width = rwidth)))
  })
  
  for (peakset in names(erl)) {
    plot_base_name <- paste0(paste(grps,collapse = "vs"),"_",peakset)
    print("plot base name:")
    print(plot_base_name)
    extended_granges <- erl[[peakset]]
    # #Keep only peaks in canonical chrs
    # canonical_chrs <- paste0("chr",c(seq(1:22),"X","Y"))
    # print("peaks in non-canonical chromosomes")
    # print(table(!seqnames(extended_granges)%in%canonical_chrs))
    # print("removing peaks in following non-canonical chrs:")
    # print(unique(as.character(seqnames(extended_granges)[!seqnames(extended_granges)%in%canonical_chrs])))
    # extended_granges <- extended_granges[seqnames(extended_granges)%in%canonical_chrs]
    # 
    ###Remove possible overlapping regions 
    tfo=findOverlaps(extended_granges,extended_granges)
    overlaps=data.frame(q=queryHits(tfo),s=subjectHits(tfo))
    overlaps <- overlaps[overlaps$q!=overlaps$s,]
    #first remove the repeated rows, keep only those in which q is smaller than s 
    overlapsf <- overlaps[overlaps$q<overlaps$s,]
    rows_to_remove <- unique(overlapsf$s)
    if(length(rows_to_remove)>0)extended_granges <- extended_granges[-rows_to_remove]
    
    
    ######Read and process coverage data####
    ##Read in coverage data and generate the variables of interest

    list_all_PeaksposCPM <- list()
    i=0
    for (bw in bwfiles) {
      #for (bw in sample_table$file[1:10]) { #to test
      i=i+1
      tmp_tssCPM=import(bw,which=extended_granges)
      #extend the GRanges object so each position around -Xkb - +X kb of each TSS has a CPM value
      all_TSSposCPM <- rep(tmp_tssCPM$score,width(ranges(tmp_tssCPM)))
      list_all_PeaksposCPM[[i]] <- all_TSSposCPM
    }
    print("Finished data import\nChecking if incomplete peaks")
    #####Check that each feature have all positions. If incomplete, remove!####
    # generate some indexes that will be useful to remove incomplete 
    #peaks in case they happen:
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
    #5) use rle to check if any incomplete peakIDs
    rle_peakIDCovobj <- rle(peakIDCovobj)
    if(any(rle_peakIDCovobj$lengths!=rwidth)){
      peaks_to_remove <- rle_peakIDCovobj$values[rle_peakIDCovobj$lengths!=rwidth]
      print("peaks to remove because incomplete:")
      print(peaks_to_remove)
      print("These correspond to")
      positions_to_remove <- peakIDCovobj%in%peaks_to_remove
      table(positions_to_remove)
      print(as.numeric(table(positions_to_remove)[2]))
      print("positions")
      #remove positions to remove in each sample position-coverage vector 
      list_all_PeaksposCPM <- lapply(list_all_PeaksposCPM,function(x)return(x[!positions_to_remove]))
      print("After removing incomplete peaks, there are left")
      print(lapply(list_all_PeaksposCPM, function(x)length(x)/rwidth))
      print("TSSs to plot")
      
    }else{
      print("There are")
      print(lapply(list_all_PeaksposCPM, function(x)length(x)/rwidth))
      print("TSSs to plot")
    }
    #####
    
    ####Compute mean CPM for each position*gene for each sample:
    mean_posCPM <- lapply(list_all_PeaksposCPM, function(x)aggregate(x,by=list(pos=rep(1:rwidth,length(x)/rwidth)),mean))
    mean_posCPM <- as.data.frame(lapply(mean_posCPM,function(s)return(s$x)))
    
    
    ######################AVERAGE PLOTS######################
    #####1) Plot density for each sample showing each group in a different color
    #define color palette (color-blind friendly)
    #define function to plot
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
    
    #test:
    #
   plot_densities(mean_posCPM,sampleTable = data.frame(group=grp),groupCol = "group",main=paste(grps,collapse = "vs"))
    
   #A) plot all samples individually
   graphics.off()
   pdf(paste0("average_plots/",plot_base_name,".pdf"))
   plot_densities(mean_posCPM,
                  sampleTable = data.frame(group=grp),
                  groupCol = "group",
                  main=paste(grps,collapse = " vs "))
   dev.off()
   
   #B) compute average for each group
   mean_posCPM_grpAvg <- t(aggregate(t(mean_posCPM),
                                            by=list(grp=grp),
                                            mean)[,-1])
   graphics.off()
   pdf(paste0("average_plots/Avg_samples_",plot_base_name,".pdf"))
   plot_densities(mean_posCPM_grpAvg,
                  sampleTable = data.frame(group=sort(grps)),
                  groupCol = "group",ylab = "mean CPM (samples avg)",
                  main=paste(grps,collapse = " vs "))
   dev.off()
   
   ###########################HEATMAPS########################
   ###2)Compute mean CPM for each position*peak for each sample set:
   #Define sample sets: has to be a factor of two groups:
   print("Computing mean CPM for each position-peak for each group")
   
   #we use it to generate average vectors for each position*peak
   #across all samples from each group
   
   list_Peak_pos_meanCPM <- lapply(grps, function(x){
     if(sum(grp==x)==1) return(as.numeric(unlist(list_all_PeaksposCPM[grp==x])))else{
       return(apply(as.data.frame(list_all_PeaksposCPM[grp==x]),1,mean))
     } 
   })
     
     # apply(as.data.frame(list_all_PeaksposCPM[sample_table[1:20,group_feat]==x]),1,mean)#to test
     
  
   
   rm(list_all_PeaksposCPM)
   gc()
   ###Convert the average vectors just created into matrices of genes vs position 
   #for this we use the functions split "split divides the data in the vector x into the groups defined by f."
   data_to_plot <- lapply(list_Peak_pos_meanCPM,function(x){
     t(as.data.frame(split(x, 
                           ceiling(seq_along(x)/rwidth))))})
   
   
   ##Next we want to transform our data in 2 different ways taking into
   #account that both matrices will be plotted together, so they need to have
   #the same gene order and the same color scale:
   #   1) sort rows in both matrices by mean total CPM signal between Cohesin and non-Cohesin
   #   2) convert values into log2 scale to plot
   #   3) generate color palette
   
   #   1) sort data to plot by mean total CPM
   total_CPMs <- lapply(data_to_plot, function(x)apply(x, 1,sum)) #for each element of our list,
   # for each peak (1 indicates rows) takes the sum of all CPMs (total CPM per row)
   #now, for each peak take the average of total CPMs per row between both groups
   mean_totalCPM <- apply(as.data.frame(total_CPMs),1,mean)
   #Finally, we use this vector to sort the rows in both matrices in the same way
   rowOrder=order(mean_totalCPM,decreasing = T)# 
   data_to_plot <- lapply(data_to_plot, function(x){
     return(x[rowOrder,])
     
   })
   
   #2)   Convert values into log2 scale to plot
   data_to_plot <- lapply(data_to_plot,function(x)return(log2(x+0.01)))
   
   plot(density(as.matrix(as.data.frame(data_to_plot))))
   summary(as.vector(as.matrix(as.data.frame(data_to_plot))))
   range_col_values=c(-6,2)
   cada=0.02
   mi=range_col_values[1]
   ma=range_col_values[2]
   # 
   
   palette.breaks=seq(mi,ma,cada) #palette breaks for changing colors
   
   palette.breaks.tmp=seq(round(min(as.data.frame(data_to_plot))-1),
                          round(max(as.data.frame(data_to_plot))+1),0.02) #palette breaks for all values
   palette.breaks.tmp <-  round(palette.breaks.tmp,2)
   color.palette.tmp=vector(length = (length(palette.breaks.tmp)-1))
   
   #I will try 3 types of palette:
   palette_list <- list(white_blue_orange=c("white","#0072B2","#D55E00"),
                        blue_white_orange=c("#0072B2","white","#D55E00"),
                        white_blue=c("white","#0072B2"))
   
   
   #Now that we have our data ready to plot we need to adjust some plotting parameters:
   
   ##Generate column labels = -1.5, 1, 0.5 , 0 (TSS), 0.5, 1, 1.5
   labcol <- rep("",rwidth)
   names(labcol)=as.character((llimit:ulimit)/1000)
   labcol[1+c(0,1,2,3,4)*round(rwidth/4)] <- names(labcol[1+c(0,1,2,3,4)*round(rwidth/4)])
   #points2label=c("-1","-0.5","0","0.5","1")
   #labcol[points2label] <- points2label
   names(data_to_plot) <- grps
   
   for (pname in names(palette_list)) {
     color.palette  <- colorRampPalette(palette_list[[pname]])(length(palette.breaks) - 1)
     #Set the palette and palette breaks so all plots have the same range of colors
     color.palette.tmp[(which(palette.breaks.tmp==as.character(mi)):which(palette.breaks.tmp==as.character(ma - cada)))]=color.palette
     color.palette.tmp[1:(which(palette.breaks.tmp==as.character(mi))-1)]=color.palette.tmp[which(palette.breaks.tmp==as.character(mi))]
     color.palette.tmp[which(palette.breaks.tmp==as.character(ma)):length(color.palette.tmp)]=color.palette.tmp[which(palette.breaks.tmp==as.character(ma -cada))]
     
     heats <- lapply(1:2,function(n){
       if(n==1) showleg=T else showleg=F
       m=data_to_plot[[n]]
       return(ComplexHeatmap::pheatmap(m,cluster_rows = F,
                                       cluster_cols = F,
                                       show_colnames = T,
                                       labels_col = labcol,
                                       show_rownames = F,
                                       color = color.palette.tmp,
                                       breaks = palette.breaks.tmp,
                                       use_raster=F,
                                       legend = showleg,
                                       fontsize_number = 9, column_names_side = c("top"),
                                       angle_col = c("0"),
                                       heatmap_legend_param = list(title="log2(CPM+0.01)"),
                                       main = names(data_to_plot)[n]))
     } )
     
     #graphics.off()
     
     # pdf(paste0(plot_base_name,".colorscale_",pname,".pdf"))
     # print(heats[[1]] + heats[[2]])
     # dev.off()
     #
     graphics.off()
     tiff(paste0("heatmaps/",plot_base_name,".colorscale_",pname,".tiff"), width = 8, height = 8, units = 'in', res = 300)
     print(heats[[1]] + heats[[2]])
     dev.off()
     
     
   }
   
  }
  }
