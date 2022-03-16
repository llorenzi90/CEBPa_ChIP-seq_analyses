#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --time=03:05:30
#SBATCH --ntasks=10

PPN=10

bamfile=$1 #example /scratch/llorenzi/ChIP_seq_Carini/MK_K27_R2/MK_K27_R2.sorted.bam

# bigwig
echo `date  +'%r'`
echo -e "################################\n\tbamcoverage started\n################################\n"

module load conda/current
conda activate deeptoolsenv

bamCoverage --numberOfProcessors $PPN --binSize 10 --normalizeUsing CPM --minMappingQuality 30 --bam $bamfile -o $bamfile.CPM.bw
echo -e "################################\n\tbamcoverage done\n################################\n"
echo `date  +'%r'`



