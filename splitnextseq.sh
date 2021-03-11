#!/bin/bash
# Julie Shay
# March 10, 2021
# This script will take a nextseq run folder with >24 samples, and split it into run folders
# with 24 samples each for input into COWBAT pipeline

while getopts 'i:o:sf:h' flag; do
	case $flag in
		i)
			INDIR=$OPTARG
			;;
		o)
			OUTPREFIX=$OPTARG
			;;
		s)
			SYMLINKS=1
			;;
		f)
			FDIR=$OPTARG
			;;
		h)
			h=1
	esac
done

checkandcopy () {
	line=$1
	tobecopied=$(ls $INDIR/$FDIR/$(echo $line | awk -F "," '{print $1}')_*.fastq.gz)
	if [ "$?" == "0" ]; then
		echo $line >> $wd/SampleSheet.csv
		for k in $tobecopied; do
			# note that a sequence name being a subset of another could cause problems here
			$copy $k $wd
		done
	else
		echo "Could not find fastq files for $(echo $line | awk -F "," '{print $1}')!"
	fi
}

spf=24 # samples per file

if [ -n "$h" ] || [ -z "$INDIR" ]; then
	echo "Usage: $0 -i INDIR -o OUTPREFIX [-s] [-f READLOC]"
	echo "This script splits a run directory into run directories with $spf samples each."
	echo "The optional -s tag will create symlinks instead of hard copies of files"
	echo "The optional -r tag allows you to specify a subdirectory of INDIR containing fastq files"
	exit
fi

if [ -n "$SYMLINKS" ]; then
	copy="ln -s"
	INDIR=$(realpath $INDIR)
else
	copy="cp -r"
fi

if [ ! -n "$FDIR" ]; then
	FDIR=""
fi

# make sure input directory exists
if [ ! -d "$INDIR" ]; then
	echo "Specified input directory does not exist!"
	exit
fi

# make sure there's a sample sheet
ss="$INDIR/SampleSheet.csv"
if [ ! -f "$ss" ]; then
	echo "Specified input directory does not have a sample sheet!"
	exit
fi

hline=$(grep -nP "^Sample_ID," $ss | sed 's/\:.*//')
scount=$(($(wc -l $ss | sed 's/\s.*//') - $hline))
dircount=$(($(($scount + $(($spf - 1))))/ $spf))

# make sure output directories don't already exist, then make output directories and copy metadata files
for (( i=1; i<=$dircount; i++)); do
	wd="${OUTPREFIX}_$i"
	if [ -d "$wd" ]; then
		echo "Specified output directory $wd already exists!"
		exit
	fi
	mkdir $wd
	for j in $(find $INDIR -maxdepth 1 -not -type d | grep -v SampleSheet.csv | grep -v .fastq.gz); do
		$copy $INDIR/$(basename $j) $wd
	done
	$copy $INDIR/InterOp $wd
	head -n $hline $ss > $wd/SampleSheet.csv
done

# now fill in the sample sheets and copy the appropriate fastq files
for (( i=1; $i<$dircount; i++)); do
	wd="${OUTPREFIX}_$i"
	for j in $(head -n $(($hline + $(($i * $spf)))) $ss | tail -n $spf); do
		checkandcopy $j
	done
done

# doing the last folder separately in case it has less than 24 samples
wd="${OUTPREFIX}_$dircount"
for j in $(tail -n $(($scount - $(($(($dircount - 1)) * $spf)))) $ss); do
	checkandcopy $j
done

echo "Finished splitting $INDIR into $dircount run folders!"
