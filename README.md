# MouseOrganogenesis_snATACseq_2020
### Code accompanying publication Pijuan-Sala et al., Nat. Cell Biol., 2020.

>  Author: Blanca Pijuan-Sala, Gottgens Lab,  University of Cambridge, UK

>  Date: February 2020

>  Description: README file explaining the scripts and the order they were executed.



**Notes:**

  - Download all the scripts and check them. You should change path to directories to your own. I have tried adding "/path/to/directory" to all those places (Sorry if I missed any!)
  - Download and install https://github.com/r3fang/snATAC
  - Place snATAC_bmat_BPS inside the bin directory of ./snATAC/bin/
  - Install:
    bedtools
    Samtools



**Directory**

wd=/path/to/directory/sample_pooled_preprocess_revision1



**Instructions**

1) Reverse complement HiSeq2500 sequencing reads to match NextSeq500

    Script: 01_HiSeq2500_revComplement.sh

    Depends on: ./bin/reverseComplement_barcodeheader_fastq.pl



2) Map the all fastq files

    We sequenced this dataset in two rounds. First, we merged the fastq for the first round:

    Script: 02_merge_Fastq_files.sh

    Then, we mapped those reads.

    Scripts:

      ⋅⋅* 03_bowtie2_map_all_together.sh

      ⋅⋅* 03_bowtie2_map_large_together.sh

      ⋅⋅* 03_bowtie2_map_small_together.sh

    Mapping scripts depend on: ./bin/map_ATACseq_v002.sh

    For the second round, we mapped it separately:
      Script: 03a_bowtie2_map_separate.sh

      Depends on:
        ./bin/map_ATACseq_v002.sh

      Output: $wd/01_SAM/new_embryo_revision1.sam


3) Merge the new SAM files with the old ones using samtools merge

    Output: $wd/01_SAM/combined.sam


4) Sort the combined SAM file
    ```
    samtools sort -O sam -T combined.sorted -o combined.sorted.sam combined.sam
    ```


    Output: $wd/01_SAM/combined.sorted.sam


5) Convert the combined & sorted SAM file to BAM
    ```
    samtools view -b -o new.sorted.bam combined.sorted.sam
```
    Output: $wd/01_BAM_BED/new.sorted.bam


6) Process into BED

    Script: 04_preprocess_all_together.sh
    Depends on:
      ./bin/preProcess_ATACseq_v004.sh

    Output folder: $wd/02_BAM_BED/

7) Filter out weird mapped reads.
```
    grep -e "^[0-9]" -e "^X" -e "^Y" embryo_revision1_nuclearGenes.bed >embryo_revision1_nuclearGenes.bed.filtered
```
    Output: $wd/02_BAM_BED/

8) Run peak calling in sample with no QC

    Script: 05_callPeaks_pooledSample_noQC.sh

    Depends on:
      ./bin/callPeaks_ATACseq_v003_ext150.sh

    Output folder: $wd/03_MACS2/


9) Run Barcode statistics

    Script: 06_barcodeStats.sh

    Depends on:
      ./bin/barcodeStats_noPeaks.sh #Please check script, as you will have to modify the path to the "genome" variable.

    Output folder: $wd/04_barcodeStats/


10) Perform nuclei QC. Here we have set promoter coverage > 0.03 and number of reads in peaks > 2000.

    Script: 07_nuclei_QC.R


11) Call peaks on sample without barcodes failing QC

    Script: 08_callPeaks_pooledSample_after_first_QC.sh
    Output: $wd/05_MACS2_afterQC/

    SLURM job:
    ```
    sbatch -o 08_callPeaks_firstQC.log --mail-user=bp382@cam.ac.uk --mail-type=END 08_callPeaks_pooledSample_after_first_QC.sh
```

12) Split peaks file into 4 files:
```
    cd $wd/05_MACS2_afterQC
    file=embryo_revision1_nuclearGenes_passedQC_summits_extended500bp_overlapped_noblacklist_TSS_merged.bed
    split -l 75855  -d $file $file'_subset'
```
13) Make matrix of barcodes that passed QC using peaks called in sample that passed QC

    Scripts:
    10_Sample_after_first_QC_makeBinMatrix_subset00.sh
    10_Sample_after_first_QC_makeBinMatrix_subset01.sh
    10_Sample_after_first_QC_makeBinMatrix_subset02.sh
    10_Sample_after_first_QC_makeBinMatrix_subset03.sh

    Depend on:
      https://github.com/r3fang/snATAC
      snATAC_bmat_BPS

    SLURM jobs:
```
    sbatch -o 10_makebinMat_00.log --mail-user=bp382@cam.ac.uk --mail-type=END 10_Sample_after_first_QC_makeBinMatrix_subset00.sh
    sbatch -o 10_makebinMat_01.log --mail-user=bp382@cam.ac.uk --mail-type=END 10_Sample_after_first_QC_makeBinMatrix_subset01.sh
    sbatch -o 10_makebinMat_02.log --mail-user=bp382@cam.ac.uk --mail-type=END 10_Sample_after_first_QC_makeBinMatrix_subset02.sh
    sbatch -o 10_makebinMat_03.log --mail-user=bp382@cam.ac.uk --mail-type=END 10_Sample_after_first_QC_makeBinMatrix_subset03.sh
```
    Output: $wd/06_matrix

14) Convert matrix into mtx binary

    Script: 11_Sample_after_first_QC_makeMTX.sh

    SLURM job:
```
    sbatch -o 11_binMat_mtx.log --mail-user=bp382@cam.ac.uk --mail-type=END 11_Sample_after_first_QC_makeMTX.sh
```
    Output: $wd/06_matrix

    Depends on:
        11_Sample_after_first_QC_makeMTX_subset00_01.R
        11_Sample_after_first_QC_makeMTX_subset02_03.R
        11_Sample_after_first_QC_makeMTX.R

15) Convert matrix into mtx raw

    Script:   11b_Sample_after_first_QC_makeMTX_raw.sh

    SLURM job:
```
    sbatch -o 11b_raw_binMat_mtx.log --mail-user=bp382@cam.ac.uk --mail-type=END 11b_Sample_after_first_QC_makeMTX_raw.sh

```
    Output: $wd/06_matrix

    Depends on:
    11b_Sample_after_first_QC_makeMTX_raw.R
    11b_Sample_after_first_QC_makeMTX_subset00_01_raw.R
    11b_Sample_after_first_QC_makeMTX_subset02_03_raw.R


16) Get peak Names:
```
    cd $wd/05_MACS2_afterQC

    sed 's/\t/_/g' embryo_revision1_nuclearGenes_passedQC_summits_extended500bp_overlapped_noblacklist_TSS_merged.bed >../06_matrix/embryo_revision1_allPeaks_passedQC_peakNames.txt
```

17) Get barcode names:
```
    cp embryo_revision1.xgi ../06_matrix/embryo_revision1_allPeaks_passedQC_barcodeNames.xgi
```
18) Make file with chromosomes not wanted (this will be used to remove them):
```
    grep -v -e "chr[0-9]" -e "chrX" -e "chrY" embryo_revision1_allPeaks_passedQC_peakNames.txt  > chr_not_wanted.txt
```
19) Generate metadata file

    Script: 12_metadata_creation_Sample_after_first_QC.R



20) Run cisTopic on binary matrix
    Script: 13b_cisTopic.sh

    Depends on:
      13b_cisTopic.R

21) Compute doublet scores on binary matrix using python script

    Script:
      14_doublet_scoring.ipynb
      HTML: 14_doublet_scoring.html


22) Remove doublets, visualise selected cells using the cistopic coordinates calculated before doublet removal and cluster nuclei:

    Script: 16_Visualisation_clustering_after_doublet_removal.ipynb
    HTML: 16_Visualisation_clustering_after_doublet_removal.html

23) Call peaks in clusters
    a. Prepare the files where each file contains the barcode names for a specific cluster.
        Script: 17_clusters_to_bigwig_preparationBefore.R


    b. Split the files containing the barcode names to allow parallelisation:
```
    for i in *.txt; do
      split -l 370  -d $i $i'_subset'
    done

    rename 's/.txt//' *

    mkdir batch_1
    mkdir batch_2
    mkdir batch_3

    ####Move the files into three folders: batch_1, batch_2, batch_3
```
    c. Grab reads containing those barcodes from the BED file
    Scripts:
      18_clusters_to_bigwig_grab_barcodes_batch_1.sh
      18_clusters_to_bigwig_grab_barcodes_batch_2.sh
      18_clusters_to_bigwig_grab_barcodes_batch_3.sh


    d. Make all the files containing the barcodes per cluster have the ".txt" at the end:
```
        for i in cluster*; do mv $i $i.txt; done
```
    e. Run script to pool files and convert to BW.

        Script: 19_clusters_to_bigwig_poolBatches_convert_toBW.sh

        SLURM job:
        ```
          sbatch -o 19_pool_BW.log --mail-user=bp382@cam.ac.uk --mail-type=END 19_clusters_to_bigwig_poolBatches_convert_toBW.sh
```
    f. Run script to call peaks:
          Script: 20_clusters_callPeaks.sh

          SLURM job:
          ```
            sbatch -o 20_callPeaks.log --mail-user=bp382@cam.ac.uk --mail-type=END 20_clusters_callPeaks.sh
```
    g. Merge peaks:
      Script: 21_merging_peaks.sh
      SLURM job:
      ```
        sbatch -o 21_mergePeaks.log --mail-user=bp382@cam.ac.uk --mail-type=END 21_merging_peaks.sh
```

24) Make matrix

    a. Split peaks file into 4 files to make computation easier:
```
    cd $wd/10_peaks_all
    file=peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist.bed
    split -l 76297  -d $file $file'_subset'
```
    b. Make matrix of barcodes that passed QC using peak set obtained above.

    Scripts:
      22_Sample_after_clusterPeakCall_makeBinMatrix_subset00.sh
      22_Sample_after_clusterPeakCall_makeBinMatrix_subset01.sh
      22_Sample_after_clusterPeakCall_makeBinMatrix_subset02.sh
      22_Sample_after_clusterPeakCall_makeBinMatrix_subset03.sh

    Depend on:
      https://github.com/r3fang/snATAC
      ./bin/snATAC_bmat_BPS

    SLURM jobs:
    ```
      sbatch -o 22_makebinMat_00.log --mail-user=bp382@cam.ac.uk --mail-type=END 22_Sample_after_clusterPeakCall_makeBinMatrix_subset00.sh
      sbatch -o 22_makebinMat_01.log --mail-user=bp382@cam.ac.uk --mail-type=END 22_Sample_after_clusterPeakCall_makeBinMatrix_subset01.sh
      sbatch -o 22_makebinMat_02.log --mail-user=bp382@cam.ac.uk --mail-type=END 22_Sample_after_clusterPeakCall_makeBinMatrix_subset02.sh
      sbatch -o 22_makebinMat_03.log --mail-user=bp382@cam.ac.uk --mail-type=END 22_Sample_after_clusterPeakCall_makeBinMatrix_subset03.sh
```

    Output: $wd/11_matrix_afterClusterQC

    c. Convert matrix into mtx binary

      Script: 23_Sample_after_clusterPeakCall_makeMTX.sh

      SLURM job:
```
      sbatch -o 23_binMat_mtx.log --mail-user=bp382@cam.ac.uk --mail-type=END 23_Sample_after_clusterPeakCall_makeMTX.sh
```
    Output: $wd/11_matrix_afterClusterQC

    Depends on:
        23_Sample_after_clusterPeakCall_makeMTX_subset00_01.R
        23_Sample_after_clusterPeakCall_makeMTX_subset02_03.R
        23_Sample_after_clusterPeakCall_makeMTX.R

    d. Convert matrix into mtx raw

      Script:   23b_Sample_after_clusterPeakCall_makeMTX_raw.sh

      SLURM job:
      ```
        sbatch -o 23b_raw_binMat_mtx.log --mail-user=bp382@cam.ac.uk --mail-type=END 23b_Sample_after_clusterPeakCall_makeMTX_raw.sh
```

        Output: $wd/11_matrix_afterClusterQC

        Depends on:
            23b_Sample_after_clusterPeakCall_makeMTX_subset00_01_raw.R
            23b_Sample_after_clusterPeakCall_makeMTX_subset02_03_raw.R
            23b_Sample_after_clusterPeakCall_makeMTX_raw.R


25) Perform Cell QC to filter based on reads in peaks

    a. Compute the number of reads in peaks.
      Script: 24_barcodeStats_readsInPeaks.sh

      SLURM job:
      ```
        sbatch -o 24_barcodeStats.log --mail-user=bp382@cam.ac.uk --mail-type=END 24_barcodeStats_readsInPeaks.sh
```
    b. Perform QC:
      Script: 25_reads_in_peaks_QC.R
      This script also adds information on the gates where nuclei were sorted.




26) Compute cisTopic on the binary matrix to visualise the data
    Script: 26b_cistopic_24.sh

    Depends on:
        26b_cistopic_readsinpeaks_24.R

    SLURM job:
    ```
    sbatch -o log_26b_cistopic_24.log --mail-user=bp382@cam.ac.uk --mail-type=END 26b_cistopic_24.sh
```

27) Get peak names and barcode names
```
      cd /path/to/directory/sample_pooled_preprocess_revision1/10_peaks_all
      sed 's/\t/_/g' peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist.bed >../11_matrix_afterClusterQC/embryo_revision1_allPeaks_passedQC_peakNames.txt

      cd /path/to/directory/sample_pooled_preprocess_revision1/04_barcodeStats

      cp embryo_revision1.xgi ../11_matrix_afterClusterQC/embryo_revision1_allPeaks_passedQC_barcodeNames.xgi
```
28) Visualise the data using cisTopic and compute clusters
    Script: 26_Final_visualisation_and_clustering.ipynb



29) Find called peaks at TSS:
```
    bedTSS=/path/to/directory/bin/Mus_musculus.GRCm38.92_TSS_minus500bp.bed
    bedFile=/path/to/directory/sample_pooled_preprocess_revision1/10_peaks_all/peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist.bed

    /path/to/directory/bin/bedtools2/bin/bedtools intersect -a $bedFile -b $bedTSS -wa -wb >peaks_celltypes_inTSS.bed
    awk '{print $1"_"$2"_"$3;}' peaks_celltypes_inTSS.bed >peaks_celltypes_inTSS.txt

```

30) Generate final metadata and plots after filtering peaks
    Script: 27_metaData_and_plots_afterPeakFilter.R
    This includes heatmap of TSS.


31) Generate bigWig tracks for each cell type
    OUTPUT: $wd/17_BigWig_celltypes/

    a. Prepare files containing cell barcodes per cell type
        Script: 29a_celltypes_to_bigwig_preparationBefore.R

    b. Split files to allow parallelisation:
```
      cd /path/to/directory/sample_pooled_preprocess_revision1/17_BigWig_celltypes/barcodes
      for i in *.txt; do
        split -l 320  -d $i $i'_subset'
      done

      rename 's/.txt//' *
      mkdir batch_1
      mkdir batch_2
      mkdir batch_3

      ####Move the files into three folders: batch_1, batch_2, batch_3
```
      c. Grab barcodes from BED file.
      Scripts:
        29b_celltypes_to_bigwig_batch_1.sh
        29b_celltypes_to_bigwig_batch_2.sh
        29b_celltypes_to_bigwig_batch_3.sh

      SLURM jobs:
      ```
        sbatch -o log_29b_batch1.log --mail-user=bp382@cam.ac.uk --mail-type=END 29b_celltypes_to_bigwig_batch_1.sh
        sbatch -o log_29b_batch2.log --mail-user=bp382@cam.ac.uk --mail-type=END 29b_celltypes_to_bigwig_batch_2.sh
        sbatch -o log_29b_batch3.log --mail-user=bp382@cam.ac.uk --mail-type=END 29b_celltypes_to_bigwig_batch_3.sh
```
      d. Make all the files containing the barcodes have the "txt" at the end:
```
          for i in cluster*; do mv $i $i.txt; done
```
      e. Convert BED files to BW:
        Script: 29c_celltypes_to_bigwig_poolBatches_convert_toBW.sh

        SLURM job:
        ```
          sbatch -o log_29c_BW.log --mail-user=bp382@cam.ac.uk --mail-type=END 29c_celltypes_to_bigwig_poolBatches_convert_toBW.sh
```
32) chromVAR on genomic regions
      Output= $wd/16_chromVAR

      a. Compute chromVAR:
          Script: 33_chromVar_computation_motifsv2_celltype_all.R

      b. Run wilcoxon rank sum tests to find cell type specific peaks:

          Script: 33b_TF_enrichment_analysis_wilcoxon_allpeaks.ipynb

      c. Get top TFs enriched.
          Script: 33b_chromVAR_celltype_specific_TFs_allpeaks.R


      d. Assess expression levels of TFs enriched in scRNA-seq dataset (Pijuan-Sala et al., Nature, 2019)
          Script: 33c_scRNAseq_E825_analysis_TFexpression.R
          This script runs along 33b_chromVAR_celltype_specific_TFs_allpeaks.R; once the top TFs are computed, this script (33c) should be run to obtain their expression levels. Then, the analysis continues on 33b, where the expression is compared with the enrichment score.

      e. Output sequence logos for TF motifs used
          Script: 33d_sequenceLogos_chromVAR.R


33) Peak annotation

    a. Generate file with no "chr"
```
        sed 's/^chr//g' peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist.bed >peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_nochr.bed
```
    b. Execute the code here to annotate regions (intergenic, intronic, exonic, TSS...)
        Script: 32b_peakAnnotation.sh
        Depends on: HOMER

    c. Match peaks to genes ensembl

      - First, take the information from Ensembl regarding gene IDs. Extract chr, start, end, strand and ensembl ID ==> mm10_geneEnsembl.txt
      https://www.ensembl.org/biomart/martview/9f07ba72335711044352afa5af3443c1

      - Run the following code:
```
        awk '{print $4"\t"$3"\t"$5"\t"$6"\t"$1"\t"$2;}' mm10_geneEnsembl.txt >mm10_geneEnsembl_sorted0.txt
        grep -e "^[0-9]" -e "^Y" -e "^X" -e "^MT" mm10_geneEnsembl_sorted0.txt >mm10_geneEnsembl_sorted.txt

        #Open it and remove first line with nano

        awk '{if($4=="1"){$3=$2;$2=$2-1000;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6;}else if($4=="-1"){$2=$3;$3=$3+1000;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6;}}' mm10_geneEnsembl_sorted.txt > mm10_geneEnsembl_sorted_gene1kbup.tmp
        awk '{if($2<1){$2=1;} if($3<1){$3=1;} print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6;}' mm10_geneEnsembl_sorted_gene1kbup.tmp > mm10_geneEnsembl_sorted_gene1kbup.txt

        #Run peaktogenes:

        cd /path/to/directory/sample_pooled_preprocess_revision1/10_peaks_all
        peak_in=peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_nochr.bed
        output=peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_nochr_geneAssignment_distance.bed
        input=mm10
        source=E
        prom_db=T

        peaktogenes.pl $peak_in $output $input $source $prom_db
        #peaktogenes.pl and its dependencies are in the bin directory. You should check the files to modify directories.


        awk '{$4=$5=$6=$7="Unmapped"; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;}' UNMAPPED_PEAKS.txt >UNMAPPED_PEAKS_mod.txt

        cat peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_nochr_geneAssignment_distance.bed UNMAPPED_PEAKS_mod.txt >peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_nochr_geneAssignment_distance_plusUnmapped.bed

```
34) Define cell type-specific regions
    Output: $wd/14_visualisation


    a. Visualise genomic regions with cells as variables
      Script: 31_Visualisation_genomic_regions_nucleiAsVars_TFIDF_in_mat.ipynb

    b. Define cell type-specific regions:
      32_Define_celltype_specific_regions.R

      One part is run in the server and another locally.
      This script also conveys all the peak annotation.
      Here you will also find the analysis of GATA binding sites in topics

      SLURM job:
```
      sbatch -o log_32_celltypepeaks.log --mail-user=bp382@cam.ac.uk --mail-type=END 32_Define_celltype_specific_regions.sh
```

35) Find regions unique to each cisTopic and perform TF enrichment analyses
    Script: 26e_cisTopic_TFmotif_unique.R

36) Check peaks of enhancers in embryoid bodies
    Scripts: 37_Heptad_peaks_overlap_bash.sh
             38_Heptad_peaks_overlap.R


***_ALLANTOIC-HAEMATO-ENDOTHELIAL LANDSCAPE_***

37) Enrichment scores for TAL1 ChIP-seq of haemangioblasts, haemogenic endothelium and HP
    a. Obtain ChIP-seq data
```
      wd=/path/to/directory/sample_pooled_preprocess_revision1/19_peaks_other_datasets
      cd $wd
      cp /path/to/directory/Experiments/PhD_BPS53/preprocess/sample_pooled_preprocess/round4/21_heptad_peaks/heptad_mm9_mm10_mod.bed ./
      mv heptad_mm9_mm10_mod.bed heptad_peaks_mm10.bed

      wget http://codex.stemcells.cam.ac.uk/data/bb/mm10/D340004_Scl.bb
      bigBedToBed D340004_Scl.bb ./TAL1_HP_peaks_mm10.bed

      wget http://codex.stemcells.cam.ac.uk/data/bb/mm10/BG251_SLX7049_Tal1.bb
      bigBedToBed BG251_SLX7049_Tal1.bb ./TAL1_HE_peaks_mm10.bed

      wget http://codex.stemcells.cam.ac.uk/data/bb/mm10/BG207_BG295_Tal1_Mesoderm.bb
      bigBedToBed BG207_BG295_Tal1_Mesoderm.bb ./TAL1_Haemangioblast_peaks_mm10.bed

      cp /path/to/directory/ChIP-Seq/GSE59402/GSM1436364-5/Peaks/GSM1436364-5_p1e-3_400bp.bed ./
      cp /path/to/directory/ChIP-Seq/GSE59402/GSM1436367-8/Peaks/GSM1436367-8_p1e-4_400bp.bed ./

      mv GSM1436364-5_p1e-3_400bp.bed ./ETV2_EB35_pAb_peaks_mm10.bed
      mv GSM1436367-8_p1e-4_400bp.bed ./ETV2_EB35_pV5_peaks_mm10.bed
      mv *.bed ./peaks_bed/

```
    b. Recompute cisTopic on landscape with allantois, endothelium and erythroid cell types
       Script: 34b_endothelium_analysis_cistopic_ery_endo_allantois.R

       SLURM job:
       ```
       sbatch -o log_34b_ery_endo_allantois.log --mail-user=bp382@cam.ac.uk --mail-type=END 34b_endothelium_analysis_cistopic_ery_endo_allantois.sh
```
    b. Compute ChIP-seq enrichment score
       Script: 34i_endothelium_all_ery_cisTopic_ChIPseq_peaks.R

       SLURM job:
       ```
       sbatch -o log_34i_cisTopic_endothelium_ChIPseq.log --mail-user=bp382@cam.ac.uk --mail-type=END 34i_endothelium_all_ery_cisTopic_ChIPseq_peaks.sh
```
38) Compute PAGA structure of allantoic-haemato-endothelial landscape
    Script: 34a_endothelium_analysis_ery_endo_allantois_publication.ipynb
    In this script, we also compute pseudotime from endothelium to allantois.


39) Find dynamically accessible regions
    Script: 34h_endothelium_trajectory_all_ery_server_v2.R



40) Find regions bound by ETV2 in dynamically accessible regions by pattern
For Regions bound by ETV2, we intersected the peaks called for the ETV2 V5 dataset with the regions assigned to each pattern.
```
    for i in {1..12}
    do
      echo "cluster $i lines"
      bedFile='/path/to/directory/sample_pooled_preprocess_revision1/18_endothelium_analysis/data/DPT_clusters_asGut/snATACseq_endothelium_eryAl_DPT_cluster_'$i'.txt.bed'
      etv2=/home/USSR/codex-pipeline/Data/Rebecca/Blanca/sample_pooled_preprocess_revision1/19_peaks_other_datasets/peaks_bed/ETV2_EB35_pV5_peaks_mm10.bed

      cat $bedFile | wc -l
      echo "cluster $i intersected"

      /home/USSR/bp382/bin/bedtools2/bin/bedtools intersect -a $bedFile -b $etv2 -wa | wc -l
    done

    #To make the plots, the numbers are typed in the script 34h_endothelium_trajectory_all_ery_server_v2.R (bottom part)
```

41) TF motif enrichment analysis using HOMER on each pattern
    a. Run commands
    ```
    cd /path/to/directory/sample_pooled_preprocess_revision1/18_endothelium_analysis/data/DPT_clusters_asGut

    for i in snATACseq_endothelium_eryAl_DPT_cluster_*.txt; do sed 's/_/\t/g' $i >$i.bed; done
```
    b. Run motif enrichment
    Script: 34m_endothelium_trajectory_all_ery_TF_enrichment.sh

    SLURM job:
    ```
    sbatch -o log_34m_DPT.log --mail-user=bp382@cam.ac.uk --mail-type=END 34m_endothelium_trajectory_all_ery_TF_enrichment.sh
```
    c. Plot heatmap of enriched motifs
    Script: 34m_endothelium_allantois_trajectory_patterns_heatmap_motifs.R


42) Check ETS factors on scRNA-seq
    a. Get the data from Pijuan-Sala et al., 2019
      Script: 35a_endothelium_scRNAseq_getData.R

    b. Perform batch correction as in Pijuan-Sala et al., 2019
      Script: 35b_endothelium_scRNAseq_batchCorrection.R

    c. Assess structure dataset and compute DPT
       35c_endothelium_scRNAseq_dataset_structure.ipynb

    d. Check expression dynamics of ETS factors
      Script: 35c_getETS_factors_testDynamics.R
