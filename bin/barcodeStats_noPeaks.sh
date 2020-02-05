#!/bin/bash

##################################################################
# calculate 3 major barcode stats
# 1) number of reads
# 2) consecutive promoter coverage
# 3) reads in peak ratio
##################################################################

if [ $# -lt 3 ]
  then
	printf "\nUsage:\n"
	printf "\nbarcodeStats.sh BEDFile promoterFile prefix\n\n"
	printf "\tBEDFile: Bed file.\n"
  printf "\tpromoterFile: File of promoters to analyse. If set to default mm10_consecutive_promoters.bed will be used.\n"
  printf "\tprefix: Prefix to use for filenames. \n"

	exit 1
fi

echo "Setting up variables..."
wd=$PWD'/'
file=$1
promFile=$2
filename=$3

if [ $promFile = "default" ]
then
	genome='/servers/bio-shares-4/gottgens/Blanca/bioinformatic_resources/mm10_consecutive_promoters.bed'
fi


# count number of reads per barcode
awk '{print $4;}' $file | sort | uniq -c | awk '{print $2"\t"$1;}' | sort -k1,1 > $filename.reads_per_cell
# consecutive promoter coverage
/home/USSR/bp382/bin/bedtools2/bin/bedtools intersect -wa -wb -a $file -b $promFile | awk '{print $4"\t"$8;}' | sort | uniq | awk '{print $1;}' | uniq -c | awk '{print $2"\t"$1;}' | sort -k1,1 > $filename.promoter_cov

