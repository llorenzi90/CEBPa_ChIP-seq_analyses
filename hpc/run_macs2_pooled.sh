#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --time=03:05:30

wdir=$1
cd $wdir
sample_id=$2
sample=$3
sample=$(echo $sample|tr ',' ' ')
control=$4
control=$(echo $control|tr ',' ' ')

outname=$sample_id.pooled.NoModel_q_0.05

module load conda/current
conda activate macs2


echo `date  +'%r'`
echo -e "################################\n\tmacs2 started for\nsample(s):$sample\ncontrol(s):$control\n################################\n"
macs2 callpeak -t $sample -c $control -n $outname -f BAMPE -g mm --qvalue 0.05 --nomodel --keep-dup 1

#when paired-end data specifying the format 'BAMPE' is required: Please note that if the format is set as BAMPE or BEDPE, MACS2 will call its special Paired-end mode to call peaks by piling up the actual ChIPed fragments defined by both aligned ends, instead of predicting the fragment size first and extending reads. Also please note that the BEDPE only contains three columns, and is NOT the same BEDPE format used by BEDTOOLS. 

#I specify --keep-dup 1 (this is the default, i.e remove duplicates) and I use a bam file with marked duplicates, although I am pretty sure this does not change anything as MACS2 calculates duplicates by itself: "Note, if you've used samtools or picard to flag reads as 'PCR/Optical duplicate' in bit 1024, MACS2 will still read them although the reads may be decided by MACS2 as duplicate later. If you plan to rely on samtools/picard/any other tool to filter duplicates, please remove those duplicate reads and save a new alignment file then ask MACS2 to keep all by '--keep-dup all'. The default is to keep one tag at the same location. Default: 1"

#q-val is also default value

echo -e "################################\n\tmacs2 done\n################################\n"
echo `date  +'%r'`

#run:
#wdir=/scratch/llorenzi/Maria_CEBPa_ChIP-seq_samples/macs2_pooled_reps

#cat macs2_pooled_sample_sheet.txt |while read s;do sample=$(echo $s|cut -d " " -f1);echo $sample; treat=$(echo $s|cut -d " " -f2);echo $treat; input=$(echo $s|cut -d " " -f3);echo $input;sbatch -J $sample.macs2_pooled run_macs2_pooled.sh $wdir $sample $treat $input;done
