#!/bin/bash
die() {
	printf '%s\n' "$1" >&2
    exit 1
    }
	

show_help() {
cat << EOF
		Usage: ${0##*/} [options] [working_dir]...
		Run ChIP-seq basic data processing in IJC hpc;
		 from QC to alignment to generation of bigWig files for visualization 
		 
		Available options are:
			
			-h			display this help and exit
			-rs runstep		run the analysis starting from this particular step. All steps are run by default (runstep=all) options: all,qc,trim,align,sort_index,markdups,filtstats,bigWig.
			-only			if this flag is set (no additional arguments) then only the step indicated by -rs will be run (default FALSE)
			-ppn integer		number of proccesors per node (default 10)
		   
EOF
}
	
runstep=all
PPN=10
only=false

#parse optional arguments	   

while :; do
	case $1 in
		-h|-\?|--help)
			show_help    # Display a usage synopsis.
			exit
			;;
		-rs|--runstep)       # Takes an option argument; ensure it has been specified.
			re='^-'
			if [ "$2" ] &&  ! [[ $2 =~ $re ]]; then
				case $2 in 
				  "all"|"qc"|"trim"|"align"|"sort_index"|"markdups"|"filtstats"|"bigWig")
					  runstep=$2
					  echo "running step $runstep"
				  ;;
				  *)
				    die 'ERROR: non-valid runnning step
					"-rs" requires one of these non-empty option arguments : "all", "qc", "trim", "align", "sort_index", "markdups", "filtstats" or "bigWig" .'
				  ;;
				esac
	            shift 2
	        else
	            die 'ERROR: "-rs" requires one of these non-empty option arguments : "all", "qc", "trim", "align", "sort_index", "markdups", "filtstats" or "bigWig" '
			fi
			;;
		-ppn)
			if [ "$2" ]; then
				PPN=$2
				re='^[0-9]+$'
				if ! [[ $PPN =~ $re ]] ; then
				   die 'ERROR:  "-ppn" must be followed by an integer' 
				fi
	            shift 2
	        else
	            die 'ERROR: "-ppn" requires an integer as option argument'
			fi
			;;
		-only)
			only=true
			shift
			echo "Only $runstep will be run"
			;;
		-?*)
			printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
			show_help
			shift
	    	;;
	    *)  # Default case: No more options, so break out of the loop.
	    	break
	esac
done	    

echo "PPN: $PPN"
echo "only: $only"
echo "wdir: $1"

if [ $runstep == all ] && $only
then
	only=false
	echo 'Note that when -rs all then -only must by FALSE, setting -only back to false'
fi
#Now the actual analysis starts 
#NOTE: this is a preliminary version, in the future I would like to make a function out of each step for the sake of organisation

sdir=$1
sample=$(basename $sdir) 


#indexes and ref files.These are mouse files, they can be changed if working with mouse or another species
bowtie2Index="/mnt/beegfs/public/references/index/bowtie2/GRCm39_UCSC/GRCm39_UCSC"
 
cd $sdir

#1) QC
step=qc
if [ $runstep == all ] || [ $runstep == $step ]
then
	read1=$(ls *_1.fq.gz)
	read2=$(ls *_2.fq.gz)
	
	if ! [ -f "$read1" ]; then
		die "input file $read1 not found"
	fi
	read2=$(ls *_2.fq.gz)
	if ! [ -f "$read2" ]; then
		die "input file $read2 not found"
	fi
	echo "read1: $read1"
	echo "read2: $read2"
	
	echo `date`
	echo -e "################################\n\tFASTQC started\n################################\n"
	module load FastQC/0.11.9
	mkdir fastqc
	fastqc $read1 -o fastqc
	fastqc $read2 -o fastqc
	echo -e "################################\n\tFASTQC done\n################################\n"
	echo `date  +'%r'`
	if $only
	then 
		exit 0
	else 
		runstep=all
	fi
fi


#2)trim_galore
step=trim
if [ $runstep == all ] || [ $runstep == $step ]
then
	read1=$(ls *_1.fq.gz)
	if ! [ -f "$read1" ]; then
		die "input file $read1 not found"
	fi
	read2=$(ls *_2.fq.gz)
	if ! [ -f "$read2" ]; then
		die "input file $read2 not found"
	fi
	echo "read1: $read1"
	echo "read2: $read2"
	
	echo `date  +'%r'`
	echo -e "################################\n\tTrimGalore started\n################################\n"

	module load TrimGalore/0.6.6


	trim_galore --paired --length 35 --fastqc $read1 $read2
	echo -e "################################\n\tTrimGalore done\n################################\n"
	echo `date  +'%r'`
	if $only
	then 
		exit 0
	else 
		runstep=all
	fi
fi
	

#3)alignment
step=align
if [ $runstep == all ] || [ $runstep == $step ]
then
	read1_val=$(ls *val_1.fq.gz)
	if ! [ -f "$read1_val" ]; then
		die "input file $read1_val not found"
	fi
	read2_val=$(ls *val_2.fq.gz)
	if ! [ -f "$read2_val" ]; then
		die "input file $read2_val not found"
	fi
	echo "read1: $read1_val"
	echo "read2: $read2_val"
	
	echo `date  +'%r'`
	echo -e "################################\n\tBowtie2 started\n################################\n"
	module load bowtie2/2.4.4
	module load gcc-11.2.0-gcc-4.8.5-pjl2s6u
	bowtie2 --very-sensitive -X 1000 -p $PPN -x $bowtie2Index -1 $read1_val -2 $read2_val -S $sample.bt2.sam
	echo -e "################################\n\tBowtie2 done\n################################\n"
	echo `date  +'%r'`

	
	if $only
	then 
		exit 0
	else 
		runstep=all
	fi
fi
	
#4) convert to bam while sorting and index
step=sort_index
if [ $runstep == all ] || [ $runstep == $step ]
then
	input=$sample.bt2.sam
	if ! [ -f "$input" ]; then
		die "input file $input not found"
	fi
	
	echo `date  +'%r'`
	echo -e "################################\n\tsamtools sort started\n################################\n"
	module load samtools-1.12-gcc-11.2.0-n7fo7p2
	
	samtools sort -@ $PPN -O bam -o $sample.sorted.bam $input
	samtools index -@ $PPN $sample.sorted.bam
	
	echo -e "################################\n\tsamtools sort done\n################################\n"
	echo `date  +'%r'`
	
	if $only
	then 
		exit 0
	else 
		runstep=all
	fi
fi



#5) mark duplicates
step=markdups
if [ $runstep == all ] || [ $runstep == $step ]
then 
	input=$sample.sorted.bam
	if ! [ -f "$input" ]; then
		die "input file $input not found"
	fi
	
	echo `date  +'%r'`
	echo -e "################################\n\tpicard Mark duplicates started\n################################\n"
	module load picard/2.26.3	
	java -jar $PICARD_DIR/picard.jar MarkDuplicates -I $input -O $sample.sorted.markedDups.bam -M $sample.dupmatrix --CREATE_INDEX true --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT --REMOVE_DUPLICATES false --TAGGING_POLICY All
	
	echo -e "################################\n\tpicard Mark duplicates done\n################################\n"
	echo `date  +'%r'`
	
	if $only
	then 
		exit 0
	else 
		runstep=all
	fi
fi


#6) MAPQ distribution, raw and proper pairs, filtering and stats
step=filtstats
if [ $runstep == all ] || [ $runstep == $step ]
then 
	input=$sample.sorted.markedDups.bam
	if ! [ -f "$input" ]; then
		die "input file $input not found"
	fi
	
	module load samtools-1.12-gcc-11.2.0-n7fo7p2
	
	echo `date  +'%r'`
	echo -e "################################\n\tsamtools view and stats\n################################\n"
	#A) FILTERS: proper pairs and MAPQ >= 2
	
	#A.1) retain only proper pairs and reads with minimum MAPQ=2
	samtools view -@ $PPN -b -h -f2 -q2 $sample.sorted.markedDups.bam > $sample.sorted.markedDups.proper_pairs.minq2.bam

		#B) MAPQ distribution for raw reads and proper pairs 
	raw=$sample.sorted.markedDups
	samtools view $raw.bam |cut -f5 |sort -n |uniq -c > $raw.MAPQdist
	samtools view -f2 $raw.bam |cut -f5 |sort -n |uniq -c > $raw.proper_pairs.MAPQdist

	#C) STATS and fragment lenght distributions for both bam files (raw and proper pairs + qulity filter) 
	for fi in $sample.sorted.markedDups $sample.sorted.markedDups.proper_pairs.minq2
	do
		samtools stats $fi.bam > $fi.stats
		samtools view $fi.bam | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > $fi.fragment_length_count.txt
	done

	#D) index the final file
	samtools index -@ $PPN $sample.sorted.markedDups.proper_pairs.minq2.bam

	#E) remove intermediate sam and bam files (only interested in the filtered one for downstream analyses)
	sam=$sample.bt2.sam
	rm $sam $sample.sorted.bam $sample.sorted.bam.bai $sample.sorted.markedDups.bam $sample.sorted.markedDups.bai 

	echo -e "################################\n\tsamtools view and stats done\n################################\n"
	echo `date  +'%r'`
	
	if $only
	then 
		exit 0
	else 
		runstep=all
	fi
fi

#7) bigwig
step=bigWig
if [ $runstep == all ] || [ $runstep == $step ]
then
  input=$sample.sorted.markedDups.primary_chr.minq2.bam
  if ! [ -f "$input" ]; then
    die "input file $input not found"
  fi

  module load deepTools/3.5.1 

  echo `date  +'%r'`
  echo -e "################################\n\tbamcoverage started\n################################\n"
  bamCoverage --binSize 10 --normalizeUsing CPM --bam $input -o $input.CPM.bw
  echo -e "################################\n\tbamcoverage done\n################################\n"
  echo `date  +'%r'`

  if $only
  then
    exit 0
  else
    runstep=all
  fi
fi

		

