#!/bin/bash

###############################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to map reads using bowtie2.
###############################################

if [ $# -lt 5 ]
  then
	printf "\nUsage:\n"
	printf "\npreProcess_ATACseq_v001.sh R1 R2 genome_dir filename num_proc\n\n"
  printf "\tExample: preProcess_ATACseq_v001.sh path/to/directory/ paired default default\n\n"
	printf "\tR1: R1 read fastq pair (unzipped).\n"
  printf "\tR2: R2 read fastq pair (unzipped).\n"
  printf "\tgenome_dir: genome directory. If default: /servers/bio-shares-4/gottgens/Blanca/bioinformatic_resources/bowtie2/genome_index/mm10_GRCm38/mm10_GRCm38 \n"
  printf "\tfilename: e.g. embryo_small_NextSeq500_1 \n"
  printf "\tnum_proc: threads \n"


	exit 1
fi


echo "Setting up variables..."
wd=$PWD'/'

R1=$1
R2=$2
genome=$3
filename=$4
num_proc=$5

if [ $genome = "default" ]
then
	genome='/path/to/directory/bioinformatic_resources/bowtie2/genome_index/mm10_GRCm38/mm10_GRCm38'
fi

cd $wd
mkdir 01_SAM
mkdir 02_BAM_BED
mkdir 03_MACS2

echo "Running bowtie2..."


bowtie2 -p 5 -t -X 2000 --no-mixed --no-discordant -p $num_proc -x $genome -1 $R1 -2 $R2 >$wd/01_SAM/$filename.sam

echo "Splitting header and alignment"

grep "^@" $wd'/01_SAM/'$filename'.sam' >$wd/01_SAM/$filename'_header.sam'

grep -v "^@" $wd'/01_SAM/'$filename'.sam' >$wd/01_SAM/$filename'_alignment.sam'


echo "DONE!"
