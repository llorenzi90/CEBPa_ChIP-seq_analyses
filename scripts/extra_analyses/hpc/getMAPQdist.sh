#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --time=02:35:30

module load samtools/1.9

bam=$1
sample=${bam/.bam/}

#1) MAPQ distribution raw

samtools view $sample.bam |cut -f5 |sort -n |uniq -c > $sample.MAPQdist


#2) MAPQ dist proper pairs
samtools view -f2 $sample.bam |cut -f5 |sort -n |uniq -c > $sample.proper_paris.MAPQdist


#3) retain only proper pairs and mapq >=2
# samtools view -h -b -f2 -q2 $sample.bam > $sample.proper_paris.minmapq2.bam
# samtools index $sample.proper_paris.minmapq2.bam

#run
#find /scratch/llorenzi/ChIP_seq_Carini/*/*sorted.markedDups.bam| while read bam; do bn=$(basename $bam); sbatch -J $bn.MAPQdist getMAPQdist.sh $bam;done 
