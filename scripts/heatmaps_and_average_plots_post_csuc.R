cArgs=commandArgs(trailingOnly = T)
options(scipen = 999) 
outdir=""
ChName=cArgs[1]
ChName="HPC7_Fli1"
outdir=paste0("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/Overlap_DBPeaks_publicChIP-seq/p30UTvsp42UT1_noTSS_wider/",ChName)

#Extend region to plot
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


dim(common_peaksdtp)

max_cols_common_peaks=apply(all_peaks_raw$common_peaks, 2, max)
summary(max_cols_common_peaks)
which.max(max_cols_common_peaks)
max(max_cols_common_peaks)
mean_cols_common_peaks=all_peaks_mean$common_peaks
which.max(mean_cols_common_peaks)
max(mean_cols_common_peaks)
plot(llimit:ulimit,max_cols_common_peaks,type = "s")
plot(llimit:ulimit,mean_cols_common_peaks,type = "s")
points(llimit:ulimit,max_cols_common_peaks,type = "s",col="red")

sort()
data_to_plot=common_peaksdtp

print(summary(as.vector(as.matrix(as.data.frame(data_to_plot)))))
plot(density(as.vector(as.matrix(as.data.frame(data_to_plot)))))
mean_cols_common_peaks=apply(common_peaksdtp, 2, mean)
which.max(mean_cols_common_peaks)


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
CHeat=ComplexHeatmap::pheatmap(m[1:100,],cluster_rows = F,
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

print(CHeat)

graphics.off()
tiff(paste0("heatmaps/",plot_base_name,".colorscale_",pname,".tiff"), width = 8, height = 8, units = 'in', res = 300)
dev.off()

rowsums=rowSums(2^data_to_plot - 0.01)
sort(rowsums,decreasing = T)==rowsums
head(rowsums)
raw_dtp=2^data_to_plot - 0.01
summary(raw_dtp[1,])
max(raw_dtp[1,])
max_vals=apply(raw_dtp,1,max)
head(max_vals)
