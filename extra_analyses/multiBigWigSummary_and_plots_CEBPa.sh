#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --time=02:05:30

module load conda/current
conda activate deeptoolsenv
cd /scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples

files=$(ls */*bw|grep -v ER)
multiBigwigSummary bins -b $files -o multiBigWig_summary_results.npz

plotCorrelation -in multiBigWig_summary_results.npz -c spearman --skipZeros -p heatmap --plotTitle "Spearman correlation based on Genome Coverage (bigWig files)" -o CEBPa_ChIP-seq_bigWig_heatmap.pdf

plotPCA -in multiBigWig_summary_results.npz --plotTitle "PCA based on Genome Coverage (bigWig files)" -o CEBPa_ChIP-seq_bigWig_PCA.pdf

