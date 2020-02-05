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
##  DESCRIPTION: Script to calculate barcode stats for reads in peaks.
###############################################################################
wd=/path/to/directory/sample_pooled_preprocess_revision1

promFile='/path/to/directory/bioinformatic_resources/mm10_consecutive_promoters.bed'


cd $wd/12_barcodeStats_celltypePeaks

bedFile=$wd/02_BAM_BED/embryo_revision1_nuclearGenes.bed.filtered.chr

echo "Setting up variables..."
file=$bedFile
filename=embryo_celltypeSpecificPeaks

peakFile=$wd'/10_peaks_all/peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist.bed'
echo "files" >file.tmp

/path/to/directory/bin/bedtools2/bin/bedtools intersect -a $file -b $peakFile -u | awk '{print $4;}' | sort | uniq -c | awk '{print $2"\t"$1;}' | sort -k1,1 - > $filename.reads_in_peak &

cd $wd/09_celltypes_Peaks
for i in *_extended500bp_qvalThreshold_overlapped.bed
do


  file=$(basename "$i")  && \
  filename="${file%_summits_extended500bp_qvalThreshold_overlapped.bed}"  && \
  filename2="${filename:49}"  && \

  echo "\n$filename start\n" && \

  /path/to/directory/bin/bedtools2/bin/bedtools intersect -a $bedFile -b $file -u | awk '{print $4;}' | sort | uniq -c | awk '{print $2"\t"$1;}' | sort -k1,1 - > ../12_barcodeStats_celltypePeaks/$filename.reads_in_peak  &&\
  echo "$filename done\n" >>../12_barcodeStats_celltypePeaks/file.tmp &
done

cd $wd/12_barcodeStats_celltypePeaks


while (( "$( ls file.tmp | wc -l | xargs)" < 22)); do
  wait
done


echo "\n...DONE!\n"
