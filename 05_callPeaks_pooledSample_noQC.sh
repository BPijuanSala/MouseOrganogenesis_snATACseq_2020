#!/bin/bash
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=1-00:00:00
#SBATCH -e slurm.MACS2alltgt.err


###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to merge BED files and to call peaks in the pooled
##               sample without having gone through Quality Control.
##               This script depends on callPeaks_ATACseq_v003_ext150.sh
##               you need to download the mm10 black list from
##               https://sites.google.com/site/anshulkundaje/projects/blacklists
###############################################################################
wd=/path/to/directory/sample_pooled_preprocess_revision1

cd $wd/03_MACS2/


bedFile=$wd/02_BAM_BED/embryo_revision1_nuclearGenes.bed.filtered




##Call peaks
sh /path/to/directory/bin/callPeaks_ATACseq_v003_ext150.sh $bedFile embryo_revision1_together_nuclearGenes


blackList=/path/to/directory/bin/mm10.blacklist.bed
bedFile=$wd/embryo_revision1_nuclearGenes_summits_extended500bp_overlapped.bed


##Remove blacklisted regions
/home/USSR/bp382/bin/bedtools2/bin/bedtools intersect -a $bedFile -b $blackList -v >embryo_revision1_nuclearGenes_summits_extended500bp_overlapped_noblacklist.bed
