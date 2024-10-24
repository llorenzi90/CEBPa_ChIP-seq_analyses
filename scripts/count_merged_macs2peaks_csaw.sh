cd "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts"

bams=$(ls ../../data/bams/|grep -v INP |grep -v ER|grep bam$)
bams=$(echo $bams |tr ' ' ',' )

outfile="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/analyses/macs2_merged_peaks_counts/merged_peaks_EV.p30.p42_LPS.UT.counts"
bedfile="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/macs2_pooled_reps/merged_peaks/merged_peaks_EV.p30.p42_LPS.UT.sorted.bed"

cd "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/bams"
Rscript "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/csaw_regionCounts_Nsamples.R" $bams "${bedfile}" "${outfile}"



