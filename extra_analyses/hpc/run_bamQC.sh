#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --time=04:35:30

module load conda/current
conda activate atacseqqc 

bampath=$1
Rscript /home/llorenzi/jobs/ChIP-seq_November2021/bamQC.R $bampath

