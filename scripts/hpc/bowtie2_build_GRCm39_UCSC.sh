#!/bin/bash
#SBATCH -o %x%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
##SBATCH -p mem
#SBATCH --ntasks=16
#SBATCH --time=03:05:30
PPN=16

module load conda/current
conda activate cutadaptenv
bowtie2="/home/llorenzi/bowtie2-2.4.2-linux-x86_64/bowtie2"
list_of_chr=$(cat /home/llorenzi/mm39_comma_separated_chr_list.txt)
cd /scratch/llorenzi/references/mouse/UCSC_mm39_genome/
  
$bowtie2-build --threads $PPN $list_of_chr GRCm39_UCSC


