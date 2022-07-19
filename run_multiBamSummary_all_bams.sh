#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -p mem #default is std
#SBATCH --time=10:05:30
#SBATCH --cpus-per-task=10

module load conda/current
conda activate deeptoolsenv

bam_files_list=$1

bam_files_one_line=$(cat $bam_files_list| tr '\n' ' ')

wdir=$(dirname $bam_files_list)
cd $wdir

multiBamSummary bins \
 --bamfiles $bam_files_one_line \ 
 --binSize 50000 \
 --minMappingQuality 2 \
 --samFlagInclude 2 \
 --smartLabels \
 -p 10 \
 -out whole_genome_readCounts.allbams.npz --outRawCounts whole_genome_readCounts.allbams.tab

#sbatch run_multiBamSummary_all_bams.sh /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/list_bam_files.txt
