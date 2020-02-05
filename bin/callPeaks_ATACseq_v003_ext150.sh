#!/bin/bash


###############################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to call and extend peaks using MACS2
###############################################

if [ $# -lt 2 ]
  then
	printf "\nUsage:\n"
	printf "\ncallPeaks_ATACseq_v003 filename.bed filename\n\n"
  printf "\tfilename.bed: bed file \n"
  printf "\tfilename: name \n"

	exit 1
fi

echo "Setting up variables..."
wd=$PWD'/'

i=$1
filename=$2

/path/to/directory/bin/python/anaconda2/bin/macs2 callpeak -t $i -f BED -n $filename -g mm -p 0.05 --nomodel --shift 0 --extsize 150 --keep-dup all


awk -v OFS="\t" '{if($2 < 250){$2=1;}else{$2=$2-250;}$3=$3+250;$1="chr"$1;print;}' $filename'_summits.bed' | sed -e 's|\s\+|\t|g' >$filename'_summits_extended500bp.bed'
/path/to/directory/bin/bedtools2/bin/bedtools merge -i $filename'_summits_extended500bp.bed' >$filename'_summits_extended500bp_overlapped.bed'

awk -v OFS="\t" '{if($2 < 37){$2=1;}else{$2=$2-37;}$3=$3+37;$1="chr"$1;print;}' $filename'_summits.bed' | sed -e 's|\s\+|\t|g' >$filename'_summits_extended74bp.bed'
/path/to/directory/bin/bedtools2/bin/bedtools merge -i $filename'_summits_extended74bp.bed' >$filename'_summits_extended74bp_overlapped.bed'


awk -v OFS="\t" '{if($2 < 500){$2=1;}else{$2=$2-500;}$3=$3+500;$1="chr"$1;print;}' $filename'_summits.bed' | sed -e 's|\s\+|\t|g' >$filename'_summits_extended1kb.bed'
/path/to/directory/bin/bedtools2/bin/bedtools merge -i $filename'_summits_extended1kb.bed' >$filename'_summits_extended1kb_overlapped.bed'
