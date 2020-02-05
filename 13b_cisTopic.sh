#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=3-00:00:00
#SBATCH -e slurm.cistopic.err

###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Execute 13b_cisTopic.R in cluster
###############################################################################

cd /path/to/directory/sample_pooled_preprocess_revision1/07_doublet_removal
/path/to/directory/R/bin/Rscript 13b_cisTopic.R
