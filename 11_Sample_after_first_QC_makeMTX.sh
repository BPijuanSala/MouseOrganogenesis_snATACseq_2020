#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=3-00:00:00
#SBATCH -e slurm.binMat_mtx_round4.err

###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Execute 11_Sample_after_first_QC_makeMTX.R in cluster
###############################################################################


/path/to/directory/R/bin/Rscript 11_Sample_after_first_QC_makeMTX_subset00_01.R
/path/to/directory/R/bin/Rscript 11_Sample_after_first_QC_makeMTX_subset02_03.R
/path/to/directory/R/bin/Rscript 11_Sample_after_first_QC_makeMTX.R
