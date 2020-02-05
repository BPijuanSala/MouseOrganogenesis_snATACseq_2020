
###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  October 2019
##  DESCRIPTION: Perform peak annotation with HOMER
##                Depends on HOMER.
###############################################################################


cd /path/to/directory/sample_pooled_preprocess_revision1/10_peaks_all

#Create annotation file
parseGTF.pl /path/to/directory/bin/Mus_musculus.GRCm38.92.gtf ann > annotations.txt
sort -k1,1 annotations.txt >annotations_sorted.txt
assignGenomeAnnotation annotations_sorted.txt annotations_sorted.txt -prioritize annotations.final.txt > stats.txt


#annotations.final.txt is the annotation file

#Add ID to peak file

awk '{print "chr"$1"_"$2"_"$3;}' peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_nochr.bed >tmp.tmp
paste peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_nochr.bed tmp.tmp >peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_nochr_ID.bed

#run genome annotation
assignGenomeAnnotation peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_nochr_ID.bed annotations.final.txt -ann peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_annotated_assignGenomeAnnotation.bed>stats2.txt

paste peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_nochr_ID.bed peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_annotated_assignGenomeAnnotation.bed >peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_annotated_assignGenomeAnnotation_withPeaks.bed
