#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -p mem #default is std
#SBATCH --time=10:05:30


#ChIP-seq & ATAC-seq data analysis from fastq file to bigwig

#set threads
#PPN=10 #this has to match ntasks 

#command line arguments
#example fastq_path:
#November 2021 ChIP-seq files: /scratch/llorenzi/ChIP_seq_Carini/HPC-7_K27_R1
dirn=$1
sample=$(basename $dirn)
#change to fastq dir
cd $dirn
#get fastq file names
#read1=$(ls *_1.fq.gz)
#read2=$(ls *_2.fq.gz)

#indexes and ref files
bowtie2Index="/scratch/llorenzi/references/mouse/UCSC_mm39_genome/bowtie2_index/GRCm39_UCSC"

#load modules and set tool paths
bowtie2="/home/llorenzi/bowtie2-2.4.2-linux-x86_64/bowtie2"
picard="/home/llorenzi/picard.jar"
trim_galore="/home/llorenzi/TrimGalore-0.6.6/trim_galore"
module load samtools/1.9
module load conda/current
#conda activate deeptoolsenv





# #1) QC
# echo `date`
# echo -e "################################\n\tFASTQC started\n################################\n"
# module load fastqc/0.11.8
# mkdir fastqc
# fastqc $fastqn -o fastqc
# echo -e "################################\n\tFASTQC done\n################################\n"
# echo `date  +'%r'`
# 
# #2)trim_galore
# echo `date  +'%r'`
# echo -e "################################\n\tTrimGalore started\n################################\n"

conda activate cutadaptenv

# $trim_galore --paired --retain_unpaired --length 35 --fastqc $read1 $read2

#file=${read1%.fq.gz}
#echo "file $file"
read1_val=$(ls *_val_1.fq.gz)
echo "read1_val $read1_val"

read2_val=$(ls *_val_2.fq.gz)
echo "read2_val $read2_val"


# echo -e "################################\n\tTrimGalore done\n################################\n"
# echo `date  +'%r'`

#3)alignment
echo `date  +'%r'`
echo -e "################################\n\tBowtie2 started\n################################\n"
$bowtie2 --very-sensitive -X 1000 -x $bowtie2Index -1 $read1_val -2 $read2_val -S $sample.bt2.sam

echo -e "################################\n\tBowtie2 done\n################################\n"
echo `date  +'%r'`

sam=$sample.bt2.sam

#4) convert to bam while sorting

echo `date  +'%r'`
echo -e "################################\n\tsamtools sort started\n################################\n"
samtools sort -O bam -o $sample.sorted.bam $sam
echo -e "################################\n\tsamtools sort done\n################################\n"
echo `date  +'%r'`

#5) index
echo `date  +'%r'`
echo -e "################################\n\tsamtools index started\n################################\n"
samtools index $sample.sorted.bam
echo -e "################################\n\tsamtools index done\n################################\n"
echo `date  +'%r'`

#6) bigwig
echo `date  +'%r'`
echo -e "################################\n\tbamcoverage started\n################################\n"
conda activate deeptoolsenv

bamCoverage --binSize 10 --normalizeUsing CPM --minMappingQuality 30 --bam $sample.sorted.bam -o $sample.sorted.bam.CPM.bw
echo -e "################################\n\tbamcoverage done\n################################\n"
echo `date  +'%r'`



#example run:
#find /scratch/llorenzi/ChIP_seq_Carini/ -maxdepth 1 -mindepth 1| while read s; do bn=$(basename $s);echo $bn;sbatch -J $bn.fastq_to_bigWig ChIP-seq_fastq_to_bigWig_NOV2021.sh $s;done
