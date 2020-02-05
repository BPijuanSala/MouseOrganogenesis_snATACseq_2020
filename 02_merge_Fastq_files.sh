#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=12:00:00
#SBATCH -e slurm.mergeFiles.err


###############################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to merge all fastq files.
###############################################


wd=/path/to/directory

cd /path/to/directory/sample_pooled_preprocess/fastq


cat $wd/NextSeq500_1/fastq/180823_NB501692_all.demultiplexed.R1.fastq $wd/NextSeq500_2/fastq/wt_e85_embryo_all.demultiplexed.R1.fastq $wd/HiSeq2500_compl/fastq/HiSeq2500_all_compl.demultiplexed.R1.fastq >embryo_all.demultiplexed.R1.fastq
cat $wd/NextSeq500_1/fastq/180823_NB501692_all.demultiplexed.R2.fastq $wd/NextSeq500_2/fastq/wt_e85_embryo_all.demultiplexed.R2.fastq $wd/HiSeq2500_compl/fastq/HiSeq2500_all_compl.demultiplexed.R2.fastq >embryo_all.demultiplexed.R2.fastq

cat $wd/NextSeq500_1/fastq/180823_NB501692_large.demultiplexed.R1.repl1.fastq $wd/NextSeq500_2/fastq/wt_e85_embryo_large.demultiplexed.R1.fastq $wd/HiSeq2500_compl/fastq/HiSeq2500_large_compl.demultiplexed.R1.fastq >embryo_large.demultiplexed.R1.fastq
cat $wd/NextSeq500_1/fastq/180823_NB501692_large.demultiplexed.R2.repl1.fastq $wd/NextSeq500_2/fastq/wt_e85_embryo_large.demultiplexed.R2.fastq $wd/HiSeq2500_compl/fastq/HiSeq2500_large_compl.demultiplexed.R2.fastq >embryo_large.demultiplexed.R2.fastq


cat $wd/NextSeq500_1/fastq/180823_NB501692_small.demultiplexed.R1.repl1.fastq $wd/NextSeq500_2/fastq/wt_e85_embryo_small.demultiplexed.R1.fastq $wd/HiSeq2500_compl/fastq/HiSeq2500_small_compl.demultiplexed.R1.fastq >embryo_small.demultiplexed.R1.fastq
cat $wd/NextSeq500_1/fastq/180823_NB501692_small.demultiplexed.R2.repl1.fastq $wd/NextSeq500_2/fastq/wt_e85_embryo_small.demultiplexed.R2.fastq $wd/HiSeq2500_compl/fastq/HiSeq2500_small_compl.demultiplexed.R2.fastq >embryo_small.demultiplexed.R2.fastq
