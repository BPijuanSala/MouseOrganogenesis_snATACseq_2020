#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=3-00:00:00
#SBATCH -e slurm.binMat_subset03.err


###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to generate binary matrix
##               Depends on https://github.com/r3fang/snATAC
##               and snATAC_bmat_BPS
###############################################################################



source /path/to/directory/bioinformatic_resources/anaconda2/bin/activate python2


wd=/path/to/directory/sample_pooled_preprocess_revision1

cd $wd/11_matrix_afterClusterQC

#Define variable for BED file
bedFile=$wd/04_barcodeStats/embryo_revision1_nuclearGenes.bed.filtered.chr.passQC
#Define variable for file containing list of barcodes to include (one barcode per line)
Barcodes=$wd/04_barcodeStats/embryo_revision1.xgi


#Define variable for file containing list of genomic regions to count
peaks=$wd/10_peaks_all/peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist.bed_subset03


echo "creating matrix "

/path/to/directory/bioinformatic_resources/snATAC/bin/snATAC_bmat_BPS -i $bedFile -x $Barcodes -y $peaks -o embryo_revision1_allPeaks_afterclusterPeak_passedQC_subset03.mat
