module load conda/current
conda activate csaw

pref=$3
samp1=$1
samp2=$2

mergePeaks -venn $pref.venn -matrix $pref.matrix -d given $samp1 $samp2 > $pref.mergedPeaks.txt
