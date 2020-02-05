#!/bin/bash


###############################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to pre-process SAM files.
###############################################


if [ $# -lt 3 ]
  then
	printf "\nUsage:\n"
	printf "\npreProcess_ATACseq_v004.sh SAM filename num_proc\n\n"
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

SAM=$1
filename=$2
num_proc=$3

cd $wd
mkdir 02_BAM_BED
mkdir 03_MACS2


echo "To BAM..."

/path/to/directory/bin/samtools-1.9_installed/bin/samtools view -@ $num_proc -bS $SAM > $wd/02_BAM_BED/$filename'.bam'

cd $wd/02_BAM_BED/
echo "preprocessing..."
#sort BAM based on readname
/path/to/directory/bin/samtools-1.9_installed/bin/samtools sort -n -O BAM -@ $num_proc $wd/02_BAM_BED/$filename'.bam' >$wd/02_BAM_BED/$filename'_sorted.bam'

/path/to/directory/bioinformatic_resources/snATAC/bin/snATAC_pre_BPS -t $num_proc -m 20 -f 2000 -e 75 -i $filename'_sorted.bam' -o $filename'.bed.gz' 2> $filename'.pre.log'
gunzip $filename'.bed.gz'
echo "Number mitochondrial reads" >>$filename"_mitochondrialReads.pre.log"
grep "MT" $filename'.bed' | wc -l >>$filename"_mitochondrialReads.pre.log"

grep -v "MT" $filename'.bed' >$filename'_nuclearGenes.bed'


echo "DONE!"
