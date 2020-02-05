#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=3-00:00:00
#SBATCH -e slurm.celltypesBW_batch3.err

###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Grab barcodes into BED files corresponding to each cluster
##               batch 3
###############################################################################

##NOTE: Before running 18_clusters_to_bigwig_grab_barcodes_batch_1.sh,
# 18_clusters_to_bigwig_grab_barcodes_batch_2.sh,
# 18_clusters_to_bigwig_grab_barcodes_batch_3.sh

# do:

#for i in *.txt; do
# split -l 320  -d $i $i'_subset'
#done
#ls *subset* | wc -l # 67
#rename 's/.txt//' *
####Move the files into three folders: batch_1, batch_2, batch_3

wd=/path/to/directory/sample_pooled_preprocess_revision1
echo "Setting up variables..."
cd $wd/08_BigWig_celltypes
rootFolder=$PWD'/'


bedFile=$wd/04_barcodeStats/embryo_revision1_nuclearGenes.bed.filtered.chr.passQC

folderBarcodes=$wd"/08_BigWig_celltypes/barcodes/batch_3/"

files=($(find $folderBarcodes -name "*subset*"))

echo "files" >file_batch3.tmp


for i in ${files[@]}; do

  file=$(basename "$i") && \
  filename=$(basename "$i") && \
  echo "\n$filename start\n" && \
  echo "\n $filename grepping barcodes\n" && \



  grep -f $i $bedFile > $bedFile"."$filename".pre" && \
  bedFile_cell=$bedFile"."$filename && \

  grep -e "^chr[0-9]" -e "^chrY" -e "^chrX" $bedFile"."$filename".pre" >$bedFile_cell && \
  echo "$filename done\n" >>file_batch3.tmp &
done



while (( "$( ls file_batch3.tmp | wc -l | xargs)" < 23)); do
  wait
done

echo "\n...DONE!\n"
