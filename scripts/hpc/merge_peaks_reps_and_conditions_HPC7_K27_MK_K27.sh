#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --time=00:35:30


cond1=HPC-7_K27
cond2=MK_K27

#generate merged peaks per replicate

cd /scratch/llorenzi/ChIP_seq_Carini
mkdir merged_peaks
cd merged_peaks

for samp in $cond1 $cond2
do
	/home/llorenzi/jobs/ChIP-seq_November2021/homer_mergePeaks_from_replicates.sh ../${samp}_R1/*narrowPeak* ../${samp}_R2/*narrowPeak* $samp.mergedReps
	mv $samp.mergedReps*R1*R2*narrowPeak $samp.mergedReps.peaks

done

/home/llorenzi/jobs/ChIP-seq_November2021/homer_mergePeaks_from_conditions.sh $cond1.mergedReps.peaks $cond2.mergedReps.peaks $cond1.$cond2
