#!/bin/bash


###############################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to calculate barcode stats: (1) Number of reads,
##               (2) constitutive promoter coverage, and (3) reads in peaks.
###############################################

if [ $# -lt 4 ]
  then
	printf "\nUsage:\n"
	printf "\nbarcodeStats.sh BEDFile promoterFile peaksFile prefix\n\n"
	printf "\tBEDFile: Bed file.\n"
  printf "\tpromoterFile: File of promoters to analyse. If set to default mm10_consecutive_promoters.bed will be used.\n"
  printf "\tpeaksFile: File with peaks to count. Normally, a narrowPeak file. BED file is also OK. \n"
  printf "\tprefix: Prefix to use for filenames. \n"

	exit 1
fi

echo "Setting up variables..."
wd=$PWD'/'
file=$1
promFile=$2
peakFile=$3
filename=$4

if [ $promFile = "default" ]
then
	genome='/path/to/directory/bioinformatic_resources/mm10_consecutive_promoters.bed'
fi


# count number of reads per barcode
awk '{print $4;}' $file | sort | uniq -c | awk '{print $2"\t"$1;}' | sort -k1,1 > $filename.reads_per_cell
# consecutive promoter coverage
/path/to/directory/bin/bedtools2/bin/bedtools intersect -wa -wb -a $file -b $promFile | awk '{print $4"\t"$8;}' | sort | uniq | awk '{print $1;}' | uniq -c | awk '{print $2"\t"$1;}' | sort -k1,1 > $filename.promoter_cov
# reads in peak ratio
/path/to/directory/bin/bedtools2/bin/bedtools intersect -a $file -b $peakFile -u | awk '{print $4;}' | sort | uniq -c | awk '{print $2"\t"$1;}' | sort -k1,1 - > $filename.reads_in_peak
