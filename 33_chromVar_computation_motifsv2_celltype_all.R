
###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to compute ChromVAR Z-scores
##               for TF motif enrichment analysis.
###############################################################################


library(SummarizedExperiment)
library(chromVAR)
library(motifmatchr)
library(BiocParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(fCI)

#######====my approach
dir <- "/path/to/directory/sample_pooled_preprocess/"

wd<- "/path/to/directory/sample_pooled_preprocess_revision1/"

source(paste0(dir,"07_motifs/scripts/chromVar_readingInputs_mod.R"))
cat("read files\n")

cat("\n reading barcodes Pass\n")

barcodesPass = as.vector(read.csv(paste0(wd,"12_barcodeStats_celltypePeaks/embryo_revision1_readsPeaks24.xgi"),
                                  sep="\n",header=F)[,1])

cat("\n reading matrix\n")

alignment_files_bed <- Matrix::readMM(paste0(wd, '11_matrix_afterClusterQC/embryo_revision1_allPeaks_afterClusterPeak_raw.mtx'))*1


cellNames = read.csv(paste0(wd,"06_matrix/embryo_revision1_allPeaks_passedQC_barcodeNames.xgi"),
                              sep="\n",header=F)

peakNames = read.csv(paste0(wd,"11_matrix_afterClusterQC/embryo_revision1_allPeaks_passedQC_peakNames.txt"),
                     sep="\n",header=F)

colnames(alignment_files_bed) = cellNames[,1]
rownames(alignment_files_bed) =peakNames[,1]

#select cells that passed QC
cellNamesSel = read.csv(paste0(wd,"12_barcodeStats_celltypePeaks/embryo_revision1_readsPeaks24.xgi"),
                     sep="\n",header=F)
alignment_files_bed = alignment_files_bed[,cellNamesSel[,1]]


#Read peaks

peakfile =  paste0(wd,'10_peaks_all/peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist.bed')
peaks0 <- getPeaks(peakfile,sort_peaks=F)
peaks <- resize(peaks0, width = 500, fix = "center")

#Read metadata

meta = read.csv(paste0(wd,"14_visualisation/data/snATACseq_embryo_revision1_metadata_afterAllQC_passQC_UMAP_cluster.txt"),
                sep="\t")
rownames(meta)=meta$barcode
meta$final_clusters = paste0("cluster_",meta$final_clusters)

cellsSelect = rownames(meta)

alignment_files_bed = alignment_files_bed[,cellsSelect]

depths = Matrix::colSums(alignment_files_bed)


cat("adding annotation clust..\n")
sample_annotation <- DataFrame(depth = depths,
                               celltype = meta$final_clusters)

counts <- SummarizedExperiment(assays = list(counts = alignment_files_bed),
                               rowRanges = peaks,
                               colData = sample_annotation)

cat("add bias\n")
counts <- chromVAR::addGCBias(counts, genome = BSgenome.Mmusculus.UCSC.mm10)

cat("add cells PassQC\n")
counts_filtered = counts[Matrix::rowSums(alignment_files_bed)>0,cellsSelect]

save(counts_filtered,file=paste0(wd,"16_chromVAR/data/CHROMVAR_embryo_all_counts_round4_peaksAll.rda"))

cat("get motifs\n")
library(chromVARmotifs)

data("mouse_pwms_v2")
motifs = mouse_pwms_v2
motif_ix <- matchMotifs(motifs, counts_filtered,
                        genome = BSgenome.Mmusculus.UCSC.mm10)

save(motif_ix,file=paste0(wd,"16_chromVAR/data/CHROMVAR_motif_ix_pwm_chromVar_peaksAll.rda"))

cat(" computing deviations\n")
dev <- computeDeviations(object = counts_filtered,
                         annotations = motif_ix)

save(dev,file=paste0(wd,"16_chromVAR/data/CHROMVAR_dev_pwm_chromVar_raw_stringentQC_refined_motifsv2_peaksAll.rda"))
devi = deviations(dev)
devi[is.na(devi)]=0

cat("saving deviations...\n")
matS <- as(devi, "sparseMatrix")
Matrix::writeMM(matS,file=paste0(wd,"16_chromVAR/data/CHROMVAR_embryo_allPeaks_jointTogether_chromVAR_deviations_sparse_pwm_chromVar_raw_stringentQC_round4_peaksAll.mtx"))
write.table(rownames(devi),file=paste0(wd,"16_chromVAR/data/CHROMVAR_embryo_allPeaks_jointTogether_chromVAR_deviations_rownames_motifs_pwm_chromVar_raw_round4_peaksAll.tab"),
            quote=F,sep="\n",col.names=F,row.names=F)
write.table(colnames(devi),file=paste0(wd,"16_chromVAR/data/CHROMVAR_embryo_allPeaks_jointTogether_chromVAR_deviations_colnames_cellNames_pwm_chromVar_raw_round4_peaksAll.tab"),
            quote=F,sep="\n",col.names=F,row.names=F)

devi = deviationScores(dev)
devi[is.na(devi)]=0

cat("saving Z-scores...\n")
matS <- as(devi, "sparseMatrix")
Matrix::writeMM(matS,file=paste0(wd,"16_chromVAR/data/CHROMVAR_embryo_allPeaks_jointTogether_chromVAR_devScores_pwm_chromVar_sparse_raw_round4_peaksAll.mtx"))
write.table(rownames(devi),file=paste0(wd,"16_chromVAR/data/CHROMVAR_embryo_allPeaks_jointTogether_chromVAR_devScores_rownames_motifs_pwm_chromVar_raw_round4_peaksAll.tab"),
            quote=F,sep="\n",col.names=F,row.names=F)
write.table(colnames(devi),file=paste0(wd,"16_chromVAR/data/CHROMVAR_embryo_allPeaks_jointTogether_chromVAR_devZscores_colnames_cellNames_pwm_chromVar_raw_round4_peaksAll.tab"),
            quote=F,sep="\n",col.names=F,row.names=F)
cat("DONE!\n\n")



motif_ix <- matchMotifs(motifs, counts_filtered,
                        genome = BSgenome.Mmusculus.UCSC.mm10,out = "positions")



save(motif_ix,file=paste0(wd,"16_chromVAR/data/CHROMVAR_motif_ix_pwm_chromVar_round4_motifv2_positioninpeak_peaksAll.rda"))
