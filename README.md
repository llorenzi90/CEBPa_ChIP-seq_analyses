# CEBPa_ChIP-seq

## The role of CEBPa in the inflammatory response 

 

The transcription factor CCAAT enhancer binding protein α (C/EBPa) has a key role in stem and progenitor cell function as it regulates the 
differentiation towards the myeloid lineage. 
As a pioneer transcription factor, it is able to bind 
nucleosome-bound DNA and associate to chromatin remodelers to 
establish new enhancers35. C/EBPa mutations are found in 15-19% of cytogenetically normal AML cases. 
The most frequent type of mutation leads to the expression of an N-terminally truncated variant of the protein, termed p30, 
instead of the full-length isoform, p42. The p30 isoform is sufficient for commitment to a myeloid fate, and the absence of p42 results in the 
transformation and expansion of the myeloid lineage. Also, it has recently been shown that loss of C/EBPa confers a competitive advantage to HSCs in 
the context of chronic inflammation. Although the consensus in the field is that acquisition of leukemic properties in C/EBPa-mutant cells is due to the 
transcriptional changes promoted by the p30 isoform, it still remains unclear what specific transcriptional programs are altered and how do they promote a 
selective advantage.  

The main mutation in C/EBPa found in AML patients consists of a frame-shift mutation that results in the aberrant expression of a shorter isoform (p30) 
instead of the normally expressed longer one (p42). To model the C/EBPa AML mutation, HPC-7 cells (which do not endogenously express C/EBPa) were transfected with 
plasmids containing either the p42 isoform or the p30 isoform. P42 and p30 are bonded to ERT2 (estrogen receptor), so that it localises in the cytoplasm until its 
induction with tamoxifen (ligand), which will translocate to the nucleus and only then p42 and p30 will be able to function. This fusion was performed to avoid 
undesired differentiative effects derived from the presence of the p42 isoform in the nucleus. We have also observed that AML patients carrying this C/EBPa
mutation have a reduced response to inflammation. To study the effects that the introduction of p30 had on HPC-7 capacity to respond to an immune challenge 
like LPS on a transcriptional level, RT-qPCR were performed. We observed that after an acute activation (2h) the expression of some acute inflammatory response 
genes like Ccl5, Cxcl2 or Il12b was not induced as much as in the cells with the long isoform of C/EBPa p42.  

## Samples  

EV 

EV + LPS 

P30 

P30 + LPS 

P42 

P42 + LPS 

 
## ChIP-seq experiments
H3k27ac ChIP-seq was performed for all samples: HPC-7 cells transfected with either p30, p42 or empty vector (EV) and subjected or not to LPS treatment (LPS vs UT).  
ER ChIP-seq was performed for p30 and p42 to check if immunoprecipitation of the bound estrogen receptor could be successfully used to pull down CEBPa.

H3K27ac is associated with the higher activation of transcription and therefore defined as an active enhancer mark. H3K27ac is found at both proximal and distal regions of transcription start site (TSS).


![image](https://user-images.githubusercontent.com/37328156/158561358-d0d66f1a-13e2-49cb-8331-d3a20e66f73e.png)



## Analyses  

### Data download 

### Run MS5 
check with file md5.txt

### QC to generation to bigWig files for visualization
ChIP-seq_fastq_to_bigWig_NOV2021_corrected_bad_indexvarname.sh

#### 1.FASTQC

#### 2.Trimming with TrimGalore!
$trim_galore --paired --retain_unpaired --length 35 --fastqc $read1 $read2

#### 3.Align to the mouse genome (mm39) with bowtie2
index building: bowtie2_build_GRCm39_UCSC.sh  
$bowtie2 --very-sensitive -X 1000 -p $PPN -x $bowtie2Index -1 $read1_val -2 $read2_val -S $sample.bt2.sam

 running Bowtie 2 with the --very-sensitive option is the same as running with options: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50. The preset options that come with Bowtie 2 are designed to cover a wide area of the speed/sensitivity/accuracy trade-off space, with the presets ending in fast generally being faster but less sensitive and less accurate, and the presets ending in sensitive generally being slower but more sensitive and more accurate. See the documentation for the preset options for details.  
 
 -D <int>
Up to <int> consecutive seed extension attempts can "fail" before Bowtie 2 moves on, using the alignments found so far. A seed extension "fails" if it does not yield a new best or a new second-best alignment. This limit is automatically adjusted up when -k or -a are specified. Default: 15.  

 -R <int>
<int> is the maximum number of times Bowtie 2 will "re-seed" reads with repetitive seeds. When "re-seeding," Bowtie 2 simply chooses a new set of reads (same length, same number of mismatches allowed) at different offsets and searches for more alignments. A read is considered to have repetitive seeds if the total number of seed hits divided by the number of seeds that aligned at least once is greater than 300. Default: 2.  
 
-N <int>
Sets the number of mismatches to allowed in a seed alignment during multiseed alignment. Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower) but increases sensitivity. Default: 0.

-L <int>
Sets the length of the seed substrings to align during multiseed alignment. Smaller values make alignment slower but more sensitive. Default: the --sensitive preset is used by default, which sets -L to 22 and 20 in --end-to-end mode and in --local mode.
 
 -i <func>
Sets a function governing the interval between seed substrings to use during multiseed alignment.  
Since it's best to use longer intervals for longer reads, this parameter sets the interval as a function of the read length, rather than a single one-size-fits-all number. For instance, specifying -i S,1,2.5 sets the interval function f to f(x) = 1 + 2.5 * sqrt(x), where x is the read length. See also: setting function options. If the function returns a result less than 1, it is rounded up to 1. Default: the --sensitive preset is used by default, which sets -i to S,1,1.15 in --end-to-end mode to -i S,1,0.75 in --local mode.  
 
 -X/--maxins <int>
The maximum fragment length for valid paired-end alignments. E.g. if -X 100 is specified and a paired-end alignment consists of two 20-bp alignments in the proper orientation with a 60-bp gap between them, that alignment is considered valid (as long as -I is also satisfied). A 61-bp gap would not be valid in that case. If trimming options -3 or -5 are also used, the -X constraint is applied with respect to the untrimmed mates, not the trimmed mates.  
The larger the difference between -I and -X, the slower Bowtie 2 will run. This is because larger differences between -I and -X require that Bowtie 2 scan a larger window to determine if a concordant alignment exists. For typical fragment length ranges (200 to 400 nucleotides), Bowtie 2 is very efficient.  
Default: 500.  
 
As of Bowtie2 v2.4.0, individual preset values can be overridden by providing the specific options e.g. the configured seed length of 20 in the [--very-senitive] preset above can be changed to 25 by also specifying the -L 25 parameter anywhere on the command line.  
 
#### 4.Convert to bam while sorting
samtools sort -@ $PPN -O bam -o $sample.sorted.bam $sam

#### 5.Samtools index
samtools index -@ $PPN $sample.sorted.bam

#### 6.Generate bigWig normalized by CPM with deeptools bamCoverage
bamCoverage --numberOfProcessors $PPN --binSize 10 --normalizeUsing CPM --minMappingQuality 30 --bam $sample.sorted.bam -o $sample.sorted.bam.CPM.bw

### Mark duplicates with Picard tools and run samtools stats
 script: picard_mark_dups.sh
 
 java -jar $picard MarkDuplicates INPUT=$sample.sorted.bam OUTPUT=$sample.sorted.markedDups.bam METRICS_FILE=$sample.dupmatrix TMP_DIR=/tmp/ CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=FALSE TAGGING_POLICY=All
 
 samtools stats $sample.sorted.markedDups.bam > $sample.sorted.markedDups.stats  
 
### Peak calling with macs2
 scripts: macs2_NoModel_q_0.05.sh and run_macs2_pooled.sh (using macs2_pooled_sample_sheet.txt)  
 Duplicates were removed in all cases (default for macs2)  
macs2 callpeak -t $sample -c $control -n $outname -f BAMPE -g mm --qvalue 0.05 --nomodel --keep-dup 1  
 
### Differential binding analysis (DESeq2) 


