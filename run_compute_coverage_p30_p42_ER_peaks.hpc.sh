#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH -p mem #default is std
##SBATCH --ntasks=10
#SBATCH --time=10:35:30
##SBATCH --mem=16G


bw_paths="/scratch/llorenzi/coverage_p30_p42/list_bigWig_files.txt"
peaks_file="/scratch/llorenzi/coverage_p30_p42/peak_sets_GR.RDS"


module load conda/current
conda activate atacseqqc

cd /scratch/llorenzi/coverage_p30_p42/

Rscript /home/llorenzi/scripts/compute_coverage_p30_p42_ER_peaks.hpc.R $bw_paths $peaks_file
