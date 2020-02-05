#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=12:00:00
#SBATCH -e slurm.prepLarge.err

###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to map fastq files from FACS gating for large size
##               nuclei using map_ATACseq_v002.sh script.
###############################################################################


wd=/path/to/directory

cd /path/to/directory/sample_pooled_preprocess/


sh /path/to/directory/bin/map_ATACseq_v002.sh $wd/sample_pooled_preprocess/fastq/embryo_large.demultiplexed.R1.fastq $wd/sample_pooled_preprocess/fastq/embryo_large.demultiplexed.R2.fastq default embryo_large 23
