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
##  DESCRIPTION: Script to pre-process aligned SAM file from indiscriminate
##               FACS gating using preProcess_ATACseq_v004.sh script.
###############################################################################

wd=/path/to/directory

#Activate Python2
source /path/to/directory/bioinformatic_resources/anaconda2/bin/activate python2

cd $wd


sh /path/to/directory/bin/preProcess_ATACseq_v004.sh $wd/sample_pooled_preprocess_revision1/01_SAM/combined.sorted.sam embryo_revision1 23
