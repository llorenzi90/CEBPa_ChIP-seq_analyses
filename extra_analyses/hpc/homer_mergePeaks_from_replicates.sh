module load conda/current
conda activate csaw

pref=$3
rep1=$1
rep2=$2

mergePeaks -venn $pref.venn -matrix $pref.matrix -d given $rep1 $rep2 -prefix $pref

