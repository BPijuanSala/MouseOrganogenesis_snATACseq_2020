#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=3-00:00:00
#SBATCH -e slurm.binMat_mtx_round4_binary.err

###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Execute 23_Sample_after_clusterPeakCall_makeMTX.R in cluster
###############################################################################


/path/to/directory/R/bin/Rscript 23_Sample_after_clusterPeakCall_makeMTX_subset00_01.R
/path/to/directory/R/bin/Rscript 23_Sample_after_clusterPeakCall_makeMTX_subset02_03.R
/path/to/directory/R/bin/Rscript 23_Sample_after_clusterPeakCall_makeMTX.R
