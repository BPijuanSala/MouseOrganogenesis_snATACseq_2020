#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=3-00:00:00
#SBATCH -e slurm.celltypesBW_peaks.err


###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: For each BED file for each cluster, call peaks
##               Depends on callPeaks_ATACseq_v003_ext150.sh
##               you need to download the mm10 black list from
##               https://sites.google.com/site/anshulkundaje/projects/blacklists
##               You also need bedtools
###############################################################################


echo "Setting up variables..."
wd=/path/to/directory/sample_pooled_preprocess_revision1

cd $wd"/09_celltypes_Peaks"
rootFolder=$PWD'/'
chrom_sizes='/path/to/directory/Programs/bioinformatics_resources/UCSC/mm10.chrom.sizes'
bedFile=$wd/04_barcodeStats/embryo_revision1_nuclearGenes.bed.filtered.chr.passQC


folderBarcodes=$wd"/08_BigWig_celltypes/barcodes/"
folderBED=$wd"/04_barcodeStats/"

files=($(find $folderBarcodes -name "*txt"))

echo "files" >file.tmp

for i in ${files[@]}; do

  file=$(basename "$i") && \
  filename="${file%.*}" && \
  dir=$(dirname "$i") && \
  echo "\n$filename start\n" && \
  #filesBED=($(find $folderBED -name "*$filename*subset0[0-9]"))
  cd $wd"/09_celltypes_Peaks/" && \
  bedFile_cell=$folderBED$filename".bed" && \
  sh /home/USSR/bp382/bin/callPeaks_ATACseq_v003_ext150.sh $bedFile_cell "embryo_revision01_nuclearGenes_passedQC_"$filename && \
  blackList=/home/USSR/bp382/bin/mm10.blacklist.bed && \
  bedFile2=$wd"/09_celltypes_Peaks/embryo_revision01_nuclearGenes_passedQC_"$filename"_summits_extended500bp_overlapped.bed" && \
  bedFile3=$wd"/09_celltypes_Peaks/embryo_all_revision01_nuclearGenes_passedQC_"$filename"_summits_extended500bp_overlapped.bed.tmp" && \
  sed 's|^chr||g' $bedFile2 > $bedFile3 && \
  mv $bedFile3 $bedFile2 && \
  /home/USSR/bp382/bin/bedtools2/bin/bedtools intersect -a $bedFile2 -b $blackList -v >"embryo_revision01_nuclearGenes_passedQC_"$filename"_summits_extended500bp_overlapped_noblacklist.bed" && \
  echo "$filename done\n" >>file.tmp &
done



while (( "$( ls file.tmp | wc -l | xargs)" < 20)); do
  wait
done

echo "\n...DONE!\n"
