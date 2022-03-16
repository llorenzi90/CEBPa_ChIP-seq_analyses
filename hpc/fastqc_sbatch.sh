#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
##SBATCH -p mem #default is std
##SBATCH --ntasks=10
##SBATCH --time=05:05:30
#SBATCH --time=00:35:30
##SBATCH --mem=16G

#ChIP-seq & ATAC-seq data analysis from fastq file to bigwig


#command line arguments
#example fastq_path:
#November 2021 ChIP-seq files: /scratch/llorenzi/ChIP_seq_Carini/HPC-7_K27_R1
dirn=$1
sample=$(basename $dirn)
#change to fastq dir
cd $dirn
#get fastq file names
read1=$(ls *_1.fq.gz|grep -v val|grep -v unpaired)
read2=$(ls *_2.fq.gz|grep -v val|grep -v unpaired)

#1) QC
echo `date`
echo -e "################################\n\tFASTQC started\n################################\n"
module load fastqc/0.11.8
mkdir fastqc
fastqc $read1 -o fastqc
fastqc $read2 -o fastqc
echo -e "################################\n\tFASTQC done\n################################\n"
echo `date  +'%r'`

#run examples:
#find /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/ -maxdepth 1 -mindepth 1 -type d|while read s; do bn=$(basename $s);echo $bn; sbatch -J $bn.fastqc fastqc.sh $s;done
#find /scratch/llorenzi/ChIP_seq_Carini/ -maxdepth 1 -mindepth 1 -type d|while read s; do bn=$(basename $s);echo $bn; sbatch -J $bn.fastqc fastqc.sh $s;done
