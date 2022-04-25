options(scipen = 999) 
region_to_plot=c(-3000,3000)
##
#lower limit
llimit=region_to_plot[1]
#upper limit
ulimit=region_to_plot[2]
#region width
rwidth=length(llimit:ulimit)


#dir.create("heatmaps")
#dir.create("average_plots")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
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

#wd="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/Overlap_DBPeaks_publicChIP-seq/p30UTvsp42UT1_noTSS_wider/"
wd="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/Overlap_DBPeaks_publicChIP-seq/p30UTvsp42UT_noTSS_wider_nogreylist/"
setwd(wd)


ChNames=list.files()
for (ChName in ChNames) {
  outdir=paste0(wd,ChName)
  print(ChName)
  setwd(outdir)
  p30_peaksdtp=readRDS(paste0(ChName,"_p30_peaks_data_to_plot_heatmap.RDS"))
  p42_peaksdtp=readRDS(paste0(ChName,"_p42_peaks_data_to_plot_heatmap.RDS"))
  common_peaksdtp=readRDS(paste0(ChName,"_common_peaks_data_to_plot_heatmap.RDS"))
  
  all_peaks=list(p30_peaks=p30_peaksdtp,
                 p42_peaks=p42_peaksdtp,
                 common_peaks=common_peaksdtp)
  
  all_peaks_raw <- lapply(all_peaks,function(x)2^x - 0.01 )
  all_peaks_mean <- lapply(all_peaks_raw,function(x)apply(x, 2, mean) )
  
  dtp=as.data.frame(all_peaks_mean)
  
  graphics.off()
  pdf(paste0("average_plots/",ChName,".pdf"))
  plot_densities(dtp,
                 sampleTable = data.frame(group=colnames(dtp)),
                 groupCol = "group",
                 cols = cbPalette[c(6,7,1)],
                 main=ChName)
  dev.off()
  
}
