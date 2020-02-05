#!/bin/bash
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=1-00:00:00
#SBATCH -e slurm.MACS2alltgt_round4.err

###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to call peaks in the pooled sample afer having
##               removed barcodes that fail Quality Control and add TSS.
##               This script depends on callPeaks_ATACseq_v003_ext150.sh
##               You need to download the mm10 black list from
##               https://sites.google.com/site/anshulkundaje/projects/blacklists
##               You also need Mus_musculus.GRCm38.92.gtf
###############################################################################


wd=/path/to/directory/sample_pooled_preprocess_revision1

cd $wd/05_MACS2_afterQC/

bedFile=$wd/02_BAM_BED/embryo_revision1_nuclearGenes.bed.filtered.chr

##Select those barcodes that have passed QC
grep -f $wd/04_barcodeStats/embryo_revision1.xgi $bedFile >$wd/04_barcodeStats/embryo_revision1_nuclearGenes.bed.filtered.chr.passQC

bedFile=$wd/04_barcodeStats/embryo_revision1_nuclearGenes.bed.filtered.chr.passQC

##Call peaks using the BED file with barcodes that have passed QC
sh /path/to/directory/bin/callPeaks_ATACseq_v003_ext150.sh $bedFile embryo_revision1_nuclearGenes_passedQC

blackList=/path/to/directory/bin/mm10.blacklist.bed
bedFile2=$wd/05_MACS2_afterQC/embryo_revision1_nuclearGenes_passedQC_summits_extended500bp_overlapped.bed
bedFile3=$wd/05_MACS2_afterQC/embryo_revision1_nuclearGenes_passedQC_summits_extended500bp_overlapped.bed.tmp

sed 's|^chr||g' $bedFile2 > $bedFile3
mv $bedFile3 $bedFile2
rm $bedFile3


##Remove blacklisted genomic regions
/home/USSR/bp382/bin/bedtools2/bin/bedtools intersect -a $bedFile2 -b $blackList -v >embryo_revision1_nuclearGenes_passedQC_summits_extended500bp_overlapped_noblacklist.bed



################
## Add TSS
################

#wd=/path/to/directory

#Make BED file from GTF file that contains TSS-500bp to TSS
#cd /path/to/directory/bin
#grep -v "^#!" Mus_musculus.GRCm38.92.gtf | awk '{if($3=="gene"){print;}}' | cut -f 9- >tmp_gtf_attributes.txt
#grep -v "^#!" Mus_musculus.GRCm38.92.gtf | awk '{if($3=="gene"){print;}}' | awk '{if($7=="+"){$5=$4;$4=$4-500;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;}else if($7=="-"){$4=$5;$5=$5+500;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;}}' > tmp_gtf_Mus_musculus.GRCm38.92_TSS_minus500bp.gtf
#grep -v "^#!" Mus_musculus.GRCm38.92.gtf | awk '{if($3=="gene"){print;}}' | awk '{print "chr"$1"\t"$4"\t"$5;}' | grep -e "chrY" -e "chrX" -e "chr[0-9]" > Mus_musculus.GRCm38.92_gene.bed
#paste tmp_gtf_Mus_musculus.GRCm38.92_TSS_minus500bp.gtf tmp_gtf_attributes.txt >tmp_gtf2_Mus_musculus.GRCm38.92_TSS_minus500bp.gtf
#grep -v "^MT" tmp_gtf2_Mus_musculus.GRCm38.92_TSS_minus500bp.gtf >Mus_musculus.GRCm38.92_TSS_minus500bp.gtf
#cut -f 9- Mus_musculus.GRCm38.92_TSS_minus500bp.gtf | awk '{print $2}' | sed -e 's/"//g' -e 's/;//g' >tmp_gtf_attributes_gene.txt
#awk '{print "chr"$1"\t"$4"\t"$5"\t";}' Mus_musculus.GRCm38.92_TSS_minus500bp.gtf >tmp_gtf_bedFile.bed
#paste tmp_gtf_bedFile.bed tmp_gtf_attributes_gene.txt >Mus_musculus.GRCm38.92_TSS_minus500bp.bed
#rm tmp_gtf*


#Make an additional file without having the gene ID.
#awk '{print $1"\t"$2"\t"$3;}' Mus_musculus.GRCm38.92_TSS_minus500bp.bed >Mus_musculus.GRCm38.92_TSS_minus500bp_noGeneName.bed


cd $wd/05_MACS2_afterQC/


bedTSS=/home/USSR/bp382/bin/Mus_musculus.GRCm38.92_TSS_minus500bp_noGeneName.bed
bedFile=$wd/05_MACS2_afterQC/embryo_revision1_nuclearGenes_passedQC_summits_extended500bp_overlapped_noblacklist.bed

#Merge the peak file with the TSS file
cat $bedFile $bedTSS | sort -k1,1 -k2,2n >embryo_revision1_nuclearGenes_passedQC_summits_extended500bp_overlapped_noblacklist_TSS.bed
/home/USSR/bp382/bin/bedtools2/bin/bedtools merge -i embryo_revision1_nuclearGenes_passedQC_summits_extended500bp_overlapped_noblacklist_TSS.bed >embryo_revision1_nuclearGenes_passedQC_summits_extended500bp_overlapped_noblacklist_TSS_merged.bed
