#!/bin/bash
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=3-00:00:00
#SBATCH -e slurm.celltypesBW_v2.err

###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Convert BED files for each cell type into BW files.
###############################################################################

echo "Setting up variables..."
wd=/path/to/directory/sample_pooled_preprocess_revision1


cd $wd"/17_BigWig_celltypes"
rootFolder=$PWD'/'

chrom_sizes='/path/to/directory/Programs/bioinformatics_resources/UCSC/mm10.chrom.sizes'
bedFile=$wd/04_barcodeStats/embryo_revision1_nuclearGenes.bed.filtered.chr.passQC

bed=$(basename "$bedFile")
folderBarcodes=$wd"/17_BigWig_celltypes/barcodes/"

folderBED=$wd"/17_BigWig_celltypes/BED/"

files=($(find $folderBarcodes -name "*txt"))

echo "\n passed\n"
echo "files" >file.tmp

for i in ${files[@]}; do

  file=$(basename "$i") && \
  filename="${file%.*}" && \
  dir=$(dirname "$i") && \
  echo "\n$filename start\n" && \
  #filesBED=($(find $folderBED -name "*$filename*subset0[0-9]"))
  cd $folderBED  && \
  echo "$filename pooling..." && \
  cat *$filename"_"*subset0[0-9] | sort -k1,1 -k2,2n >$filename".bed" && \
  cd $wd"/17_BigWig_celltypes" && \

  bdgFile=$wd"/17_BigWig_celltypes/embryo_revision1_passQCAll_nuclearGenes_"$filename".bdg" && \
  bdgFile_process=$wd"/17_BigWig_celltypes/embryo_revision1_passQCAll_nuclearGenes_processed_"$filename".bdg" && \
  bwFile_process=$wd"/17_BigWig_celltypes/embryo_revision1_passQCAll_nuclearGenes_processed_"$filename".bw" && \
  bdgFile_process_norm=$wd"/17_BigWig_celltypes/embryo_revision1_passQCAll_nuclearGenes_processed_norm_"$filename".bdg" && \
  bwFile_process_norm=$wd"/17_BigWig_celltypes/embryo_revision1_passQCAll_nuclearGenes_processed_norm_"$filename".bw" && \
  echo "\n $filename grepping barcodes\n" && \
  bedFile_cell=$folderBED$filename".bed" && \
  echo "\n$filename To Bdg...\n" && \
  /home/USSR/codex-pipeline/Programs/bedtools2/bin/bedtools genomecov -bg -i $bedFile_cell -g $chrom_sizes  >$bdgFile && \
  grep -e "^chr[0-9]" -e "^chrY" -e "^chrX" $bdgFile | sort -k1,1 -k2,2n >$bdgFile_process && \
  echo "\n$filename To Bw...\n" && \
  /home/USSR/bp382/bin/bedGraphToBigWig $bdgFile_process $chrom_sizes $bwFile_process && \
  awk 'BEGIN{SUM=0;}FNR==NR{SUM+=$4;next;}{$4=(($4/SUM)*10e7);print;}' $bdgFile_process $bdgFile_process | sed -e 's|\s\+|\t|g' | sort -k1,1 -k2,2n >$bdgFile_process_norm && \
  echo "$filename To Bw..." && \
  /home/USSR/bp382/bin/bedGraphToBigWig $bdgFile_process_norm $chrom_sizes $bwFile_process_norm  && \
  echo "$filename done\n" >>file.tmp &
done



while (( "$( ls file.tmp | wc -l | xargs)" < 20)); do
  wait
done

echo "\n...DONE!\n"
