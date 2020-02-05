#!/bin/bash
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=2-00:00:00
#SBATCH -e slurm.BCstats.err

###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to calculate barcode stats: (1) Number of reads,
##               (2) constitutive promoter coverage, and (3) reads in peaks.
##               It uses barcodeSats.sh
###############################################################################

wd=/home/USSR/codex-pipeline/Data/Blanca/sample_pooled_preprocess_revision1
promFile='/servers/bio-shares-4/gottgens/Blanca/bioinformatic_resources/mm10_consecutive_promoters.bed'

bedFile=$wd/02_BAM_BED/embryo_revision1_nuclearGenes.bed.filtered

#Add "chr" to coordinates of BED file (check if it's necessary in your case)
awk '{$1="chr"$1; print;}' $bedFile | sed -e 's|\s\+|\t|g'  >$bedFile.chr


cd $wd/04_barcodeStats

#Execute barcodeStats.sh
sh /home/USSR/bp382/bin/barcodeStats_noPeaks.sh $bedFile.chr $promFile embryo_revision1
