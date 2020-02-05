#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=3-00:00:00
#SBATCH -e slurm.binMat_subset02.err


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


cd $wd/06_matrix

#Define variable for BED file
bedFile=$wd/04_barcodeStats/embryo_revision1_nuclearGenes.bed.filtered.chr.passQC
#Define variable for file containing list of barcodes to include (one barcode per line)
Barcodes=$wd/04_barcodeStats/embryo_revision1.xgi


#Define variable for file containing list of genomic regions to count
peaks=$wd/05_MACS2_afterQC/embryo_revision1_nuclearGenes_passedQC_summits_extended500bp_overlapped_noblacklist_TSS_merged.bed_subset02

echo "creating matrix "

/path/to/directory/bioinformatic_resources/snATAC/bin/snATAC_bmat_BPS -i $bedFile -x $Barcodes -y $peaks -o embryo_revision1_allPeaks_passedQC_subset02.mat
