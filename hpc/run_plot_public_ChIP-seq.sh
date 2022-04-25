#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH -p mem #default is std
#SBATCH --ntasks=10
#SBATCH --time=5:35:30
##SBATCH --mem=16G

bw=$1
ChName=$2

module load conda/current
conda activate atacseqqc

echo $ChName
echo -e "Rscript '/home/llorenzi/scripts/plot_coverage_public_ChIP-seq.hpc.R' ${bw} ${ChName}"
Rscript /home/llorenzi/scripts/plot_coverage_public_ChIP-seq.hpc.R $bw $ChName

#run : cat ChIPs_to_plot.txt |while read bw cn; do sbatch -J $cn.plots run_plot_public_ChIP-seq.sh $bw $cn ; done


