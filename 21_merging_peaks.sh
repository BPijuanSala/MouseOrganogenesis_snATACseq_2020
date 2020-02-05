#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=3-00:00:00
#SBATCH -e slurm.celltypes_merge_peaks.err


###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Merge peaks in pooled sample + TSS + cluster-wise peaks
###############################################################################


wd=/path/to/directory/sample_pooled_preprocess_revision1

cd $wd/09_celltypes_Peaks

for i in *.narrowPeak; do

    file=$(basename "$i")
    filename="${file%_peaks.*}"
    narrowPeak=$filename"_peaks.narrowPeak"
    summits=$filename"_summits.bed"
    summitsPre=$filename"_summits"

    awk '{print $9;}' $narrowPeak >qval.tmp
    paste $summits qval.tmp >$summits.qval
    awk '{if ($6 > 30){print;}}' $summits.qval >$summits.qval.thres
    rm qval.tmp

   awk -v OFS="\t" '{if($2 < 250){$2=1;}else{$2=$2-250;}$3=$3+250;$1=$1;print;}' $summits.qval.thres | sed -e 's|\s\+|\t|g' | grep -e "chrX" -e "chrY" -e "chr[0-9]" >$summitsPre'_extended500bp_qvalThreshold.bed'
   /path/to/directory/bin/bedtools2/bin/bedtools merge -i $summitsPre'_extended500bp_qvalThreshold.bed' >$summitsPre'_extended500bp_qvalThreshold_overlapped.bed'

done

cat *_extended500bp_qvalThreshold_overlapped.bed | sort -k1,1 -k2,2n >../10_peaks_all/peaks_celltypes_extended500bp_qvalThreshold_overlapped.bed


cd $wd/10_peaks_all
cat peaks_celltypes_extended500bp_qvalThreshold_overlapped.bed ../05_MACS2_afterQC/embryo_revision1_nuclearGenes_passedQC_summits_extended500bp_overlapped_noblacklist_TSS_merged.bed | sort -k1,1 -k2,2n >peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks.bed

/path/to/directory/bin/bedtools2/bin/bedtools merge -i peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks.bed | grep -e "chrX" -e "chrY" -e "chr[0-9]" >peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged.bed


blackList=/path/to/directory/bin/mm10.blacklist.bed


bedFile=peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged.bed


/path/to/directory/bin/bedtools2/bin/bedtools intersect -a $bedFile -b $blackList -v >peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist.bed
