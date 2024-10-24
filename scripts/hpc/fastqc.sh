#!/bin/bash

module load fastqc/0.11.8

input=$1
dn=$(dirname $input)

echo -e "################################\n\tFASTQC started\n################################\n"

case $input in
    *.fq|*.fastq|*.fq.gz|*.fastq.gz)  
	mkdir $dn/fastqc
        echo "input is a fastq file"
	echo -e "\nrunning:"
	echo "fastqc $input -o $dn/fastqc"
	fastqc $input -o $dn/fastqc
        echo -e "################################\n\tFASTQC done\n################################\n"
        ;;
    *)
        echo "Input is not a fastq file, assuming it is a directory"
	mkdir $input/fastqc
	fastq_files=$(ls $input/*.{fastq,fq,fq.gz,fastq.gz})
	for fq in $fastq_files
	do
	    echo -e "\nrunning:"
	    echo "fastqc $fq -o $input/fastqc"
	    fastqc $fq -o $input/fastqc
	done
	echo -e "################################\n\tFASTQC done\n################################\n"
        ;;
esac

