###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Calculate intersections between ChIP-seq peak of
##               Tal1 (Scl) and heptad with endothelium-specific peaks.
###############################################################################


cd /path/to/directory/sample_pooled_preprocess_revision1/19_peaks_other_datasets
endo=/path/to/directory/sample_pooled_preprocess_revision1/20_celltype_specific_peaks/peaks/cluster_9.txt.bed
HE=/path/to/directory/sample_pooled_preprocess_revision1/19_peaks_other_datasets/peaks_bed/TAL1_HE_peaks_mm10.bed
heptad=/path/to/directory/sample_pooled_preprocess_revision1/19_peaks_other_datasets/peaks_bed/heptad_peaks_mm10.bed

/home/USSR/bp382/bin/bedtools2/bin/bedtools intersect -a $endo -b $HE -wa >intersected.tmp
cat intersected.tmp | wc -l
/home/USSR/bp382/bin/bedtools2/bin/bedtools intersect -a intersected.tmp  -b $heptad -wa | wc -l
