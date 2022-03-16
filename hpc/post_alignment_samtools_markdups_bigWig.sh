#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --time=3:05:30
#SBATCH --ntasks=10

PPN=${2-10} #if not passed as second argument PPN=10 

dirn=$1

sample=$(basename $dirn)
sam=$sample.bt2.sam

cd $dirn

#load modules and set tool paths

picard="/home/llorenzi/picard.jar"

module load samtools/1.9
module load conda/current
#conda activate deeptoolsenv

#4) convert to bam while sorting

echo `date  +'%r'`
echo -e "################################\n\tsamtools sort started\n################################\n"
samtools sort -@ $PPN -O bam -o $sample.sorted.bam $sam
echo -e "################################\n\tsamtools sort done\n################################\n"
echo `date  +'%r'`


#5) index
echo `date  +'%r'`
echo -e "################################\n\tsamtools index started\n################################\n"
samtools index -@ $PPN $sample.sorted.bam
echo -e "################################\n\tsamtools index done\n################################\n"
echo `date  +'%r'`

#6) mark duplicates
PATH=/home/llorenzi/jdk-16.0.1/bin:$PATH
export PATH

echo `date  +'%r'`
echo -e "################################\n\tpicard Mark duplicates started\n################################\n"
java -jar $picard MarkDuplicates INPUT=$sample.sorted.bam OUTPUT=$sample.sorted.markedDups.bam METRICS_FILE=$sample.dupmatrix TMP_DIR=/tmp/ CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=FALSE TAGGING_POLICY=All
echo -e "################################\n\tpicard Mark duplicates done\n################################\n"
echo `date  +'%r'`


#7) MAPQ distribution raw and proper pairs and stats
echo -e "################################\n\tsamtools view MAPQ distribution\n################################\n"
samtools view $sample.sorted.bam |cut -f5 |sort -n |uniq -c > $sample.sorted.MAPQdist
samtools view -f2 $sample.sorted.bam |cut -f5 |sort -n |uniq -c > $sample.sorted.proper_paris.MAPQdist
samtools stats $sample.sorted.markedDups.bam > $sample.sorted.markedDups.stats
echo -e "################################\n\tsamtools view MAPQ distribution done\n################################\n"
echo `date  +'%r'`

#8) bigwig
echo `date  +'%r'`
echo -e "################################\n\tbamcoverage started\n################################\n"
conda activate deeptoolsenv

bamCoverage --binSize 10 --normalizeUsing CPM --minMappingQuality 2 --bam $sample.sorted.bam -o $sample.sorted.bam.CPM.bw
echo -e "################################\n\tbamcoverage done\n################################\n"
echo `date  +'%r'`




