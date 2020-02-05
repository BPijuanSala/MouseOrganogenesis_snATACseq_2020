cd /home/USSR/bp382/bin


grep -v "^#!" Mus_musculus.GRCm38.92.gtf | awk '{if($3=="gene"){print;}}' | cut -f 9- >tmp_gtf_attributes.txt


grep -v "^#!" Mus_musculus.GRCm38.92.gtf | awk '{if($3=="gene"){print;}}' | awk '{if($7=="+"){$5=$5;$4=$4-1000;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;}else if($7=="-"){$4=$4;$5=$5+1000;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;}}' > tmp_gtf_Mus_musculus.GRCm38.92_gene1kbup.gtf
paste tmp_gtf_Mus_musculus.GRCm38.92_gene1kbup.gtf tmp_gtf_attributes.txt >tmp_gtf2_Mus_musculus.GRCm38.92_gene1kbup.gtf

grep -v "^MT" tmp_gtf2_Mus_musculus.GRCm38.92_gene1kbup.gtf >Mus_musculus.GRCm38.92_gene1kbup.gtf

cut -f 9- Mus_musculus.GRCm38.92_gene1kbup.gtf | awk '{print $2}' | sed -e 's/"//g' -e 's/;//g' >tmp_gtf_attributes_gene.txt

awk '{print "chr"$1"\t"$4"\t"$5"\t";}' Mus_musculus.GRCm38.92_gene1kbup.gtf >tmp_gtf_bedFile.bed
paste tmp_gtf_bedFile.bed tmp_gtf_attributes_gene.txt >Mus_musculus.GRCm38.92_gene1kbup_geneName.bed


rm tmp_gtf*



############################

cd /home/USSR/bp382/bin


grep -v "^#!" Mus_musculus.GRCm38.92.gtf | awk '{if($3=="gene"){print;}}' | cut -f 9- >tmp_gtf_attributes.txt


grep -v "^#!" Mus_musculus.GRCm38.92.gtf | awk '{if($3=="gene"){print;}}' | awk '{if($7=="+"){$5=$4+250;$4=$4-250;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;}else if($7=="-"){$4=$5-250;$5=$5+250;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;}}' > tmp_gtf_Mus_musculus.GRCm38.92_TSS_plusminus250bp.gtf
paste tmp_gtf_Mus_musculus.GRCm38.92_TSS_plusminus250bp.gtf tmp_gtf_attributes.txt >tmp_gtf2_Mus_musculus.GRCm38.92_TSS_plusminus250bp.gtf

grep -v "^MT" tmp_gtf2_Mus_musculus.GRCm38.92_TSS_plusminus250bp.gtf >Mus_musculus.GRCm38.92_TSS_plusminus250bp.gtf

cut -f 9- Mus_musculus.GRCm38.92_TSS_plusminus250bp.gtf | awk '{print $2}' | sed -e 's/"//g' -e 's/;//g' >tmp_gtf_attributes_gene.txt

awk '{print "chr"$1"\t"$4"\t"$5"\t";}' Mus_musculus.GRCm38.92_TSS_plusminus250bp.gtf >tmp_gtf_bedFile.bed
paste tmp_gtf_bedFile.bed tmp_gtf_attributes_gene.txt >Mus_musculus.GRCm38.92_TSS_plusminus250bp.bed


rm tmp_gtf*




#########


cd /home/USSR/bp382/bin


grep -v "^#!" Mus_musculus.GRCm38.92.gtf | awk '{if($3=="gene"){print;}}' | cut -f 9- >tmp_gtf_attributes.txt


grep -v "^#!" Mus_musculus.GRCm38.92.gtf | awk '{if($3=="gene"){print;}}' | awk '{if($7=="+"){$5=$4;$4=$4-500;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;}else if($7=="-"){$4=$5;$5=$5+500;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;}}' > tmp_gtf_Mus_musculus.GRCm38.92_TSS_minus500bp.gtf
paste tmp_gtf_Mus_musculus.GRCm38.92_TSS_minus500bp.gtf tmp_gtf_attributes.txt >tmp_gtf2_Mus_musculus.GRCm38.92_TSS_minus500bp.gtf

grep -v "^MT" tmp_gtf2_Mus_musculus.GRCm38.92_TSS_minus500bp.gtf >Mus_musculus.GRCm38.92_TSS_minus500bp.gtf

cut -f 9- Mus_musculus.GRCm38.92_TSS_minus500bp.gtf | awk '{print $2}' | sed -e 's/"//g' -e 's/;//g' >tmp_gtf_attributes_gene.txt

awk '{print "chr"$1"\t"$4"\t"$5"\t";}' Mus_musculus.GRCm38.92_TSS_minus500bp.gtf >tmp_gtf_bedFile.bed
paste tmp_gtf_bedFile.bed tmp_gtf_attributes_gene.txt >Mus_musculus.GRCm38.92_TSS_minus500bp.bed


rm tmp_gtf*


###############
awk '{if($4=="1"){$5=$5;$4=$4-1000;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6;}else if($7=="-1"){$4=$4;$5=$5+1000;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6;}}' > mm10_geneEnsembl_sorted_gene1kbup.gtf
