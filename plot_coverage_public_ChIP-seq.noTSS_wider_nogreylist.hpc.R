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
deseq2file="/scratch/llorenzi/HPC7_public_ChIP_plots/p30UTvsp42UT_DESeq2_res.noNAs.annotated.csv"
GLfile="/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/DiffBind/greylisted_peaks_overlap_manualapproach.peakID.txt"
greylistpeaks=read.table(GLfile)
#list of ChIP-datasets:
#list_ChIPs=list.files("/scratch/llorenzi/HPC7_public_ChIP_plots/bigWig_files/")
# ChIPName <- gsub(".bt2|[_-]ChIPSeq|[_-]ChIP-seq|[_-]ChIP_seq","",
#      gsub(".sorted.bam.CPM.bw","",basename(list_ChIPs)))
#list_ChIPs=data.frame(list_ChIPs,ChIPName)

##NOTE: these bigWigs were generated with a bowtie2 index 
#with chromosome ids instead of chromosome names,
#so we need to translate them in order to import the 
#desired granges

chr_conversion <- read.table("/scratch/llorenzi/HPC7_public_ChIP_plots/mm39_alias.tab")
colnames(chr_conversion) <- c("sequenceName",	"alias", "names",	"UCSC database: mm39")

chr_conversion$`UCSC database: mm39`[chr_conversion$sequenceName=="chrX"] <- chr_conversion$names[chr_conversion$sequenceName=="chrX"]
chr_conversion$`UCSC database: mm39`[chr_conversion$sequenceName=="chrY"] <- chr_conversion$names[chr_conversion$sequenceName=="chrY"]


#comparison p30 vs p42 (UT)
deseq2res=read.csv(deseq2file)
anyNA(deseq2res$padj)

#Remove promoter peaks:
deseq2res <- deseq2res[!grepl("Promoter",deseq2res$annotation),]

#Remove peaks in DiffBind greylist
deseq2res <- deseq2res[!deseq2res$X%in%greylistpeaks$V1,]

#extract sets of peaks
padjco=0.05
p30_peaks=deseq2res[deseq2res$padj<=padjco&deseq2res$log2FoldChange>0,]
p42_peaks=deseq2res[deseq2res$padj<=padjco&deseq2res$log2FoldChange<0,]
common_peaks=deseq2res[deseq2res$padj>padjco,]

#Extend region to plot
region_to_plot=c(-3000,3000)
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


#parse command arguments
cArgs=commandArgs(trailingOnly = T)
bw=cArgs[1]
ChName=cArgs[2]


if(length(cArgs)==2) what_to_plot="both" else{
  if(cArgs[3]==1|cArgs[3]=="avg") what_to_plot="avg" else if(cArgs[3]==2|cArgs[3]=="hmap"){
    what_to_plot="hmap"
  } else what_to_plot="both"
  
}
if(any(cArgs=="-o")|any(cArgs=="-outdir")){
  outdir=cArgs[which(cArgs%in%c("-o","-outdir"))+1]
}else outdir="/scratch/llorenzi/HPC7_public_ChIP_plots/p30UTvsp42UT_noTSS_wider_nogreylist/"

# if(any(cArgs=="-d")|any(cArgs=="-datadir")){
#   datadir=cArgs[which(cArgs%in%c("-d","-datadir"))+1]
# }else datadir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundENCODEPeaks/RDS_files/"


######Read and process coverage data####
##Read in coverage data and generate the variables of interest

    
setwd(outdir)
dir.create(ChName)
setwd(ChName)
dir.create("heatmaps")
dir.create("average_plots")
    
list_mean_posCPM <- list()
for (peakset in names(erl)) {
  plot_base_name <- paste0(ChName,"_",peakset)
  print("plot base name:")
  print(plot_base_name)
  extended_granges <- erl[[peakset]]
  
  #if extended_ranges is too long (longer than round(2001*5000/rwidth)) then split into chunks    
  chunk_length=round(10005000/rwidth)
  chunks=list()
  if(ceiling(length(extended_granges)/chunk_length)>1){
    chunks=ceiling(length(extended_granges)/chunk_length)
    cl=ceiling(length(extended_granges)/chunks)
    
    i=1
    j=cl
    chunks=list()
    n=1
    while (j<length(extended_granges)) {
      chunks[[n]]=i:j
      i=j+1
      j=i+cl -1
      n=n+1
    }
    chunks[[n]]=i:length(extended_granges)
  }else chunks[[1]]=1:length(extended_granges)
  
  all_TSSposCPM=c()
  for (cn in 1:length(chunks)) {
    ext_gr=extended_granges[chunks[[cn]]]
    tmp_tssCPM=import(bw,which=ext_gr)
    #extend the GRanges object so each position around -Xkb - +X kb of each TSS has a CPM value
    tmp_all_TSSposCPM <- rep(tmp_tssCPM$score,width(ranges(tmp_tssCPM)))
    
    print(paste0("Finished data import chunk ",cn, " of ", length(chunks)))
    print("Checking if incomplete peaks")
    #####Check that each feature have all positions. If incomplete, remove!####
    
    posIDorigin <- paste0(as.character(rep(seqnames(ext_gr),each=rwidth)),"_",
                          unlist(apply(as.data.frame(ranges(ext_gr)),1,function(x)return(x[1]:x[2]))))
    #2)generate corresponding peakIDs for each position in the original peaks:      
    chrs_origin <- as.character(rep(seqnames(ext_gr),each=rwidth))
    ranges_origin <- rep(unlist(apply(as.data.frame(ranges(ext_gr)),
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
      tmp_all_TSSposCPM <- tmp_all_TSSposCPM[!positions_to_remove]
      print("After removing incomplete peaks, there are left")
      print(length(tmp_all_TSSposCPM)/rwidth)
      print("peaks to plot")
      
    }else{
      print("No incomplete peaks. There are")
      print(length(tmp_all_TSSposCPM)/rwidth)
      print("peaks to plot")
    }
    
    all_TSSposCPM=c(all_TSSposCPM,tmp_all_TSSposCPM)
    
  }
  
  if(what_to_plot=="both"|what_to_plot=="hmap"){
  
  ###########################HEATMAPS########################
  
  ###Convert the CPM per position in a data frame
  #for this we use the functions split "split divides the data in the vector x into the groups defined by f."
  print("Converting vector of CPM per position into CPM data frame of peaks x distance from peak center")
  data_to_plot <- t(as.data.frame(split(all_TSSposCPM,
                                        ceiling(seq_along(all_TSSposCPM)/rwidth))))
  
  
  ##Next we want to transform our data in 2 different ways taking into
  #account that all ChIPs will be plotted together, so they need to have
  #the same color scale:
  #   1) sort rows in by total CPM signal
  #   2) convert values into log2 scale to plot
  #   3) generate color palette to be used for all ChIPs
  
  #   1) sort peaks to plot by total CPM
  total_CPMs <- rowSums(data_to_plot) 
  #Finally, we use this vector to sort the rows 
  rowOrder=order(total_CPMs,decreasing = T)#
  data_to_plot <- data_to_plot[rowOrder,]
  
  #2)   Convert values into log2 scale to plot
  data_to_plot <- log2(data_to_plot+0.01)
  saveRDS(data_to_plot,paste0(plot_base_name,"_data_to_plot_heatmap.RDS"))
  #plot(density(as.matrix(as.data.frame(data_to_plot))))
  print(summary(as.vector(as.matrix(as.data.frame(data_to_plot)))))
  
  range_col_values=c(-7,2) 
  cada=0.02
  mi=range_col_values[1]
  if(mi<round(min(as.data.frame(data_to_plot))-1))mi=round(min(as.data.frame(data_to_plot))-1)
  ma=range_col_values[2]
  if(ma>round(max(as.data.frame(data_to_plot))+1))ma=round(max(as.data.frame(data_to_plot))+1)
  
  
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
  
  
  #for (pname in names(palette_list)) {
  pname="blue_white_orange"
    color.palette  <- colorRampPalette(palette_list[[pname]])(length(palette.breaks) - 1)
    #Set the palette and palette breaks so all plots have the same range of colors
    color.palette.tmp[(which(palette.breaks.tmp==as.character(mi)):which(palette.breaks.tmp==as.character(ma - cada)))]=color.palette
    color.palette.tmp[1:(which(palette.breaks.tmp==as.character(mi))-1)]=color.palette.tmp[which(palette.breaks.tmp==as.character(mi))]
    color.palette.tmp[which(palette.breaks.tmp==as.character(ma)):length(color.palette.tmp)]=color.palette.tmp[which(palette.breaks.tmp==as.character(ma -cada))]
    
    print("generating tiff plot...")
    showleg=T
    m=data_to_plot
    CHeat=ComplexHeatmap::pheatmap(m,cluster_rows = F,
                                   cluster_cols = F,
                                   show_colnames = T,
                                   labels_col = labcol,
                                   show_rownames = F,
                                   color = color.palette.tmp,
                                   breaks = palette.breaks.tmp,
                                   use_raster=T,
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
    
  #}
  
  print("Computing mean CPM across all peaks for each position relative to peak center...")
  ####Compute mean CPM for each position*gene for each sample:
  mean_posCPM <- aggregate(all_TSSposCPM,by=list(pos=rep(1:rwidth,length(all_TSSposCPM)/rwidth)),mean)
  #mean_posCPM <- as.data.frame(lapply(mean_posCPM,function(s)return(s$x)))
  list_mean_posCPM[[peakset]] <- mean_posCPM
  
  
}

if(what_to_plot=="both"|what_to_plot=="avg"){
  
saveRDS(list_mean_posCPM,paste0(ChName,"_data_to_plot_average.RDS"))

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
    
    
   