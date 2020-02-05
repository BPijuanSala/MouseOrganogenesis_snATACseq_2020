#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=3-00:00:00
#SBATCH -e slurm.celltypesBW_batch2.err


###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Grab barcodes into BED files corresponding to each cell type
##               batch 2
###############################################################################

##NOTE: Before running 29_celltypes_to_bigwig_batch_1.sh,
# 29_celltypes_to_bigwig_batch_2.sh,
# 29_celltypes_to_bigwig_batch_3.sh

# do:

#for i in *.txt; do
# split -l 320  -d $i $i'_subset'
#done
#ls *subset* | wc -l # 67
#rename 's/.txt//' *
####Move the files into three folders: batch_1, batch_2, batch_3


echo "Setting up variables..."
wd=/path/to/directory/sample_pooled_preprocess_revision1
cd $wd"/17_BigWig_celltypes/BED/"
rootFolder=$PWD'/'

chrom_sizes='/path/to/directory/Programs/bioinformatics_resources/UCSC/mm10.chrom.sizes'
bedFile=$wd/04_barcodeStats/embryo_revision1_nuclearGenes.bed.filtered.chr.passQC

bed=$(basename "$bedFile")

folderBarcodes=$wd"/17_BigWig_celltypes/barcodes/batch_2/"


files=($(find $folderBarcodes -name "*subset*"))

echo "files" >file_batch2.tmp

for i in ${files[@]}; do

  file=$(basename "$i") && \
  filename=$(basename "$i") && \
  #filename="${file%.*}" && \
  echo "\n$filename start\n" && \
  echo "\n $filename grepping barcodes\n" && \



  grep -f $i $bedFile > $bed"."$filename".pre" && \
  bedFile_cell=$bed"."$filename && \
  grep -e "^chr[0-9]" -e "^chrY" -e "^chrX" $bed"."$filename".pre" >$bedFile_cell && \
  echo "$filename done\n" >>file_batch2.tmp &
done



while (( "$( ls file_batch2.tmp | wc -l | xargs)" < 25)); do
  wait
done

echo "\n...DONE!\n"
