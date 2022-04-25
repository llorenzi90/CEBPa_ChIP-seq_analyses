## ---------------------------
##
##
## Purpose of script: plot coverage of public ChIP around 
## differential peaks
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-03-29
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
require(rtracklayer)
## ---------------------------
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Define function to plot average plots:
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

#Differential peaks folder:
deseq2file='/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/DESeq2_results/annotated_results/p30UTvsp42UT_DESeq2_res.noNAs.annotated.csv' 

#list of ChIP-datasets:
list_ChIPs=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/HPC7_public_ChIP_seq/data/bigWig_files/ChIPbigWigFiles_to_plot.txt",sep = "\t")

list_ChIPs$ChIPName <- gsub(".bt2|[_-]ChIPSeq|[_-]ChIP-seq|[_-]ChIP_seq","",
     gsub(".sorted.bam.CPM.bw","",basename(list_ChIPs$V1)))

##NOTE: these bigWigs were generated with a bowtie2 index 
#with chromosome ids instead of chromosome names,
#so we need to translate them in order to import the 
#desired granges

chr_conversion <- read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/mouse/genomes/mm39_alias.tab")
colnames(chr_conversion) <- c("sequenceName",	"alias", "names",	"UCSC database: mm39")

chr_conversion$`UCSC database: mm39`[chr_conversion$sequenceName=="chrX"] <- chr_conversion$names[chr_conversion$sequenceName=="chrX"]
chr_conversion$`UCSC database: mm39`[chr_conversion$sequenceName=="chrY"] <- chr_conversion$names[chr_conversion$sequenceName=="chrY"]


#comparison p30 vs p42 (UT)
deseq2res=read.csv(deseq2file)
anyNA(deseq2res$padj)

#extract sets of peaks
padjco=0.05
p30_peaks=deseq2res %>% filter(padj<=padjco,log2FoldChange>0)
p42_peaks=deseq2res %>% filter(padj<=padjco,log2FoldChange<0)
common_peaks=deseq2res %>% filter(padj>padjco)

region_to_plot=c(-1000,1000)
##
#lower limit
llimit=region_to_plot[1]
#upper limit
ulimit=region_to_plot[2]
#region width
rwidth=length(llimit:ulimit)

pl <- list(p30_peaks=p30_peaks,
           p42_peaks=p42_peaks,
           common_peaks=common_peaks)
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
  extended_granges <-  GRanges(chrs,
                  IRanges(start = midpoints + llimit,
                          width = rwidth))
  ###Remove possible overlapping regions 
  tfo=findOverlaps(extended_granges,extended_granges)
  overlaps=data.frame(q=queryHits(tfo),s=subjectHits(tfo))
  overlaps <- overlaps[overlaps$q!=overlaps$s,]
  #first remove the repeated rows, keep only those in which q is smaller than s 
  overlapsf <- overlaps[overlaps$q<overlaps$s,]
  rows_to_remove <- unique(overlapsf$s)
  if(length(rows_to_remove)>0)extended_granges <- extended_granges[-rows_to_remove]
  
  #translate extended granges chr names to UCSC ids
  chrnames=as.character(seqnames(extended_granges))
  trchrs=chr_conversion$`UCSC database: mm39`[match(chrnames,
                                                    chr_conversion$sequenceName)]
  
  seqlevels(extended_granges) <- levels(as.factor(trchrs))
  seqnames(extended_granges) <- trchrs
  return(extended_granges)
})



  ######Read and process coverage data####
  ##Read in coverage data and generate the variables of interest
 
 
for (ch in 1:nrow(list_ChIPs)) {
    
  setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/Overlap_DBPeaks_publicChIP-seq/p30UTvsp42UT/")
  ChName=list_ChIPs[ch,2]
  dir.create(ChName)
  setwd(ChName)
  dir.create("heatmaps")
  dir.create("average_plots")
  bw=list_ChIPs[ch,1]
    
  list_mean_posCPM <- list()
  for (peakset in names(erl)) {
    plot_base_name <- paste0(ChName,"_",peakset)
    print("plot base name:")
    print(plot_base_name)
      
    extended_granges <- erl[[peakset]]
    
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
      
    
    
    if(length(extended_granges)>5000){
      
    }
    tmp_tssCPM=import(bw,which=extended_granges)
    #extend the GRanges object so each position around -Xkb - +X kb of each TSS has a CPM value
    all_TSSposCPM <- rep(tmp_tssCPM$score,width(ranges(tmp_tssCPM)))
    
    print("Finished data import\nChecking if incomplete peaks")
    #####Check that each feature have all positions. If incomplete, remove!####
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
      all_TSSposCPM <- all_TSSposCPM[!positions_to_remove]
      print("After removing incomplete peaks, there are left")
      print(length(all_TSSposCPM)/rwidth)
      print("peaks to plot")
      
    }else{
      print("There are")
      print(length(all_TSSposCPM)/rwidth)
      print("peaks to plot")
    }
    

    ###########################HEATMAPS########################
    ###2)Compute mean CPM for each position*peak for each sample set:
    #Define sample sets: has to be a factor of two groups:
    print("Computing mean CPM for each position-peak for each group")

    #we use it to generate average vectors for each position*peak
    #across all samples from each group

    ###Convert the average vectors just created into matrices of genes vs position
    #for this we use the functions split "split divides the data in the vector x into the groups defined by f."
    data_to_plot <- t(as.data.frame(split(all_TSSposCPM,
                            ceiling(seq_along(all_TSSposCPM)/rwidth))))


    ##Next we want to transform our data in 2 different ways taking into
    #account that all ChIPs will be plotted together, so they need to have
    #the same color scale:
    #   1) sort rows in by total CPM signal
    #   2) convert values into log2 scale to plot
    #   3) generate color palette to be used for all ChIPs

    #   1) sort peaks to plot by total CPM
    total_CPMs <- rowSums(data_to_plot) #for each element of our list,
    # for each peak (1 indicates rows) takes the sum of all CPMs (total CPM per row)
    #now, for each peak take the average of total CPMs per row between both groups
    #Finally, we use this vector to sort the rows in both matrices in the same way
    rowOrder=order(total_CPMs,decreasing = T)#
    data_to_plot <- data_to_plot[rowOrder,]

    #2)   Convert values into log2 scale to plot
    data_to_plot <- log2(data_to_plot+0.01)

    plot(density(as.matrix(as.data.frame(data_to_plot))))
    summary(as.vector(as.matrix(as.data.frame(data_to_plot))))
    range_col_values=c(-7,2)
    cada=0.02
    mi=range_col_values[1]
    if(mi<round(min(as.data.frame(data_to_plot))-1))mi=round(min(as.data.frame(data_to_plot))-1)
    ma=range_col_values[2]
    if(ma>round(max(as.data.frame(data_to_plot))+1))ma=round(max(as.data.frame(data_to_plot))+1)
    #

    palette.breaks=seq(mi,ma,cada) #palette breaks for changing colors

    palette.breaks.tmp=seq(round(min(as.data.frame(data_to_plot))-1),
                           round(max(as.data.frame(data_to_plot))+1),cada) #palette breaks for all values
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


    for (pname in names(palette_list)) {
      color.palette  <- colorRampPalette(palette_list[[pname]])(length(palette.breaks) - 1)
      #Set the palette and palette breaks so all plots have the same range of colors
      color.palette.tmp[(which(palette.breaks.tmp==as.character(mi)):which(palette.breaks.tmp==as.character(ma - cada)))]=color.palette
      color.palette.tmp[1:(which(palette.breaks.tmp==as.character(mi))-1)]=color.palette.tmp[which(palette.breaks.tmp==as.character(mi))]
      color.palette.tmp[which(palette.breaks.tmp==as.character(ma)):length(color.palette.tmp)]=color.palette.tmp[which(palette.breaks.tmp==as.character(ma -cada))]


        showleg=T
        m=data_to_plot
        CHeat=ComplexHeatmap::pheatmap(m,cluster_rows = F,
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
                                        main = ChName)

        graphics.off()
        tiff(paste0("heatmaps/",plot_base_name,".colorscale_",pname,".tiff"), width = 8, height = 8, units = 'in', res = 300)
        print(CHeat)
        dev.off()

    }
    
    ####Compute mean CPM for each position*gene for each sample:
    mean_posCPM <- aggregate(all_TSSposCPM,by=list(pos=rep(1:rwidth,length(all_TSSposCPM)/rwidth)),mean)
    #mean_posCPM <- as.data.frame(lapply(mean_posCPM,function(s)return(s$x)))
    list_mean_posCPM[[peakset]] <- mean_posCPM
      
    
  }
  
  
  
  ######################AVERAGE PLOTS######################
  #####1) Plot density for each sample showing each group in a different color
  #define color palette (color-blind friendly)
  #define function to plot
  
  dtp=as.data.frame(lapply(list_mean_posCPM,function(x)return(x$x)))
  colnames(dtp) <- names(list_mean_posCPM)
  
  graphics.off()
  pdf(paste0("average_plots/",ChName,".pdf"))
  plot_densities(dtp,
                 sampleTable = data.frame(group=colnames(dtp)),
                 groupCol = "group",
                 cols = cbPalette[c(6,7,1)],
                 main=ChName)
  dev.off()
  
  
      
    }
    
    
   