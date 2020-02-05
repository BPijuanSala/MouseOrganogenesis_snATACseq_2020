#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=12:00:00
#SBATCH -e slurm.allPreprocess.err


###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to map fastq files from indiscriminate FACS gating
##               using map_ATACseq_v002.sh script.
###############################################################################


wd=/path/to/directory

cd /path/to/directory/sample_pooled_preprocess/


sh /path/to/directory/bin/map_ATACseq_v002.sh $wd/sample_pooled_preprocess/fastq/embryo_all.demultiplexed.R1.fastq $wd/sample_pooled_preprocess/fastq/embryo_all.demultiplexed.R2.fastq default embryo_all 23
