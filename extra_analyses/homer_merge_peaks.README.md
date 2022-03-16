### homer merge peaks
These analyses were done because we wanted to have a quick comparison of called peaks between conditions

From homer doc:   

**Separating Peaks into Unique and Overlapping sets**   
Merging peaks together into a single file is very useful for certain types of analysis, such as making scatter plots that compare the tag-densities between peaks from separate experiments - in this case you want to count tags at specific and common regions.  Alternatively, you may be interested in separating the peaks into common and specific sets for focused analysis.  To do this use the "-prefix <filename>" option - this will create separate files based on overlapping peaks for each set of peaks. For example:  

mergePeaks -d 100 pu1.peaks cebp.peaks -prefix mmm

This will create files named "mmm_pu1.peaks", "mmm_cebp.peaks", and "mmm_pu1.peaks_cebp.peaks".

The output file will contain the following columns:

1. Merged Peak name (will start with "Merged-")
2. chromosome
3. start (average from merged peaks)
4. end (average from merged peaks)
5. strand
6. Average peak score (actually, the average of the original values in column 6 of the peak files - or column 5 of BED files)
7. Original peak files contributing to the merged peak
8. Total number of peaks merged (occasionally more than one peak from a single file will be merged if the peaks are within the specify distance or two or more peaks from one file overlap with the same single peak(s) from another file)


  **Peak Co-Occurrence Statistics**  
The mergePeaks program will also find calculate the statistics of co-occurrence between peaks in a pairwise fashion.  If "-matrix <filename>" is specified, HOMER will calculate statistics about the pairwise overlap of peaks.  Three separate pairwise matrix files will be produced using the supplied <filename> as a prefix:

filename.logPvalue.matrix.txt (natural log p-values for overlap using the hypergeometric distribution, positive values signify divergence)
filename.logRatio.matrix.txt (natural log of the ratio of observed overlapping peaks to the expected number of overlapping peaks)
filename.count.matrix.txt (raw counts of overlapping peaks)
 
The statistics are dependent on the effective size of the genome, which can be specified using "-gsize <#>" (default: 2,000,000,000)
