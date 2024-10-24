#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --time=03:05:30

sample=$1
control=$2

outname=$sample.NoModel_q_0.05

module load conda/current
conda activate macs2

echo `date  +'%r'`
echo -e "################################\n\tmacs2 started for\nsample:$sample\ncontrol:$control\n################################\n"
macs2 callpeak -t $sample -c $control -n $outname -f BAMPE -g mm --qvalue 0.05 --nomodel --keep-dup 1

#when paired-end data specifying the format 'BAMPE' is required: Please note that if the format is set as BAMPE or BEDPE, MACS2 will call its special Paired-end mode to call peaks by piling up the actual ChIPed fragments defined by both aligned ends, instead of predicting the fragment size first and extending reads. Also please note that the BEDPE only contains three columns, and is NOT the same BEDPE format used by BEDTOOLS. 

#I specify --keep-dup 1 (this is the default, i.e remove duplicates) and I use a bam file with marked duplicates, although I am pretty sure this does not change anything as MACS2 calculates duplicates by itself: "Note, if you've used samtools or picard to flag reads as 'PCR/Optical duplicate' in bit 1024, MACS2 will still read them although the reads may be decided by MACS2 as duplicate later. If you plan to rely on samtools/picard/any other tool to filter duplicates, please remove those duplicate reads and save a new alignment file then ask MACS2 to keep all by '--keep-dup all'. The default is to keep one tag at the same location. Default: 1"

#q-val is also default value

echo -e "################################\n\tmacs2 done\n################################\n"
echo `date  +'%r'`

#run:
#cat sample_control_macs2_Carini.txt |while read s; do myarray=($s); sample=${myarray[0]};echo "sample $sample"; control=${myarray[1]}; echo -e "control $control";bn=$(basename $sample);bn=${bn/.sorted.markedDups.bam/};echo $bn;sbatch -J $bn.macs2 macs2_NoModel_q_0.05.sh $sample $control ;done
#cat sample_control_macs2_Carini.txt:
#/scratch/llorenzi/ChIP_seq_Carini/MK_K27_R2/MK_K27_R2.sorted.markedDups.bam /scratch/llorenzi/ChIP_seq_Carini/INPUT_MK_R2/INPUT_MK_R2.sorted.markedDups.bam
#/scratch/llorenzi/ChIP_seq_Carini/HPC-7_K27_R1/HPC-7_K27_R1.sorted.markedDups.bam /scratch/llorenzi/ChIP_seq_Carini/INPUT_HPC-7_R1/INPUT_HPC-7_R1.sorted.markedDups.bam
#/scratch/llorenzi/ChIP_seq_Carini/MK_GATA1_R2/MK_GATA1_R2.sorted.markedDups.bam /scratch/llorenzi/ChIP_seq_Carini/INPUT_MK_R2/INPUT_MK_R2.sorted.markedDups.bam
#/scratch/llorenzi/ChIP_seq_Carini/MK_K27_R1/MK_K27_R1.sorted.markedDups.bam /scratch/llorenzi/ChIP_seq_Carini/INPUT_MK_R1/INPUT_MK_R1.sorted.markedDups.bam
#/scratch/llorenzi/ChIP_seq_Carini/HPC-7_K27_R2/HPC-7_K27_R2.sorted.markedDups.bam /scratch/llorenzi/ChIP_seq_Carini/INPUT_HPC-7_R2/INPUT_HPC-7_R2.sorted.markedDups.bam

