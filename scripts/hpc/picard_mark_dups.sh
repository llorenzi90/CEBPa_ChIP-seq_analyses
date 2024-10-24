#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --time=03:05:30

#November 2021 ChIP-seq files: /scratch/llorenzi/ChIP_seq_Carini/HPC-7_K27_R1
dirn=$1
sample=$(basename $dirn)
#change to sample dir
cd $dirn

# Mark duplicates with picard
echo `date  +'%r'`
echo -e "### $sample picard Mark duplicates ###" 

PATH=/home/llorenzi/jdk-16.0.1/bin:$PATH
export PATH
picard="/home/llorenzi/picard.jar"

java -jar $picard MarkDuplicates INPUT=$sample.sorted.bam OUTPUT=$sample.sorted.markedDups.bam METRICS_FILE=$sample.dupmatrix TMP_DIR=/tmp/ CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=FALSE TAGGING_POLICY=All

#run samtools stats on marked dups
module load samtools/1.9
samtools stats $sample.sorted.markedDups.bam > $sample.sorted.markedDups.stats

echo -e "### $sample picard Mark duplicates and related stats done###" 
echo `date  +'%r'`




