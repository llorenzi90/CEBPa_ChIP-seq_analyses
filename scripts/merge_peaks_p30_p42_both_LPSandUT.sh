#!/bin/bash

cd "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/macs2_pooled_reps"

samp=$1 #p30 o p42

cut -f1-4,9,6 ${samp}_UT.pooled.NoModel_q_0.05_peaks.narrowPeak > merged_peaks/${samp}_LPS.UT_q_0.05_peaks.narrowPeak 
cut -f1-4,9,6 ${samp}_LPS.pooled.NoModel_q_0.05_peaks.narrowPeak >> merged_peaks/${samp}_LPS.UT_q_0.05_peaks.narrowPeak 

cd merged_peaks

sort -k1,1 -k2,2n ${samp}_LPS.UT_q_0.05_peaks.narrowPeak > ${samp}_LPS.UT_q_0.05_peaks.sorted.narrowPeak

rm ${samp}_LPS.UT_q_0.05_peaks.narrowPeak

sed "s/.pooled.NoModel_q_0.05_peak_[0-9a-z]*//" ${samp}_LPS.UT_q_0.05_peaks.sorted.narrowPeak | bedtools merge -c 4,4,6 -o count_distinct,distinct,mean -i - > merged_peaks.${samp}_LPS.UT_q_0.05_peaks.sorted.narrowPeak.bed

#25/02/2022 
#Now I add also the merging of all peaks together to have a common peak file to perform quantification and differential binding across all samples 
cd /home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/CEBPa/ChIP-seq/data/peaks/macs2/macs2_pooled_reps
ls |grep narrowPeak| grep "LPS\\|UT" > list_peaks_files_to_merge.txt #peaks to merge: macs2 called on pooled reps for:
#EV_LPS
#EV_UT
#p30_LPS
#p30_UT
#p42_LPS
#p42_UT

#Merging:
files2merge=$(ls|grep narrowPeak| grep "LPS\\|UT")
cat $files2merge |cut -f1-4,9,6 | sort -k1,1 -k2,2n | sed "s/.pooled.NoModel_q_0.05_peak_[0-9a-z]*//"| bedtools merge -c 4,4,6 -o count_distinct,distinct,mean -i - > merged_peaks_EV.p30.p42_LPS.UT.sorted.bed


