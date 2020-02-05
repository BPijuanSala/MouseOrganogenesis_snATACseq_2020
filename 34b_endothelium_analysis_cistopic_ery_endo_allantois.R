###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  October 2019
##  DESCRIPTION: Recompute cisTopic on allantois, endothelium and erythroid cell types
###############################################################################


suppressWarnings(library(cisTopic))
library(Matrix)

wd = "/path/to/directory/sample_pooled_preprocess_revision1/"

cat("\n reading matrix\n")
m <- Matrix::readMM(paste0(wd, '11_matrix_afterClusterQC/embryo_revision1_allPeaks_afterClusterPeak.mat.bin.mtx'))*1

cat(dim(m))

cat("\n reading barcodes\n")

barcodes = as.vector(read.csv(paste0(wd,"04_barcodeStats/embryo_revision1.xgi"),
                              sep="\n",header=F)[,1])   


cat(dim(barcodes))
colnames(m) <- barcodes
cat("\n reading bed file\n")

regions_frame = read.table(paste0(wd,"10_peaks_all/peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist.bed"),
                           header=F, sep="\t")
colnames(regions_frame) <- c('Chr', 'Start', 'End')

rownames(m) <- paste0(regions_frame[,1], ':', regions_frame[,2], '-', regions_frame[,3])

cat("\n reading barcodes Pass\n")
umapScanpy<-read.csv(paste0(wd,"14_visualisation/data/snATACseq_embryo_revision1_reads_in_peaks_above_24_UMAP.csv"))
clustScanpy<-read.csv(paste0(wd,"14_visualisation/data/snATACseq_embryo_revision1_reads_in_peaks_above_24_Louvain.csv"))
rownames(umapScanpy)<-clustScanpy[,1]

meta = data.frame(
  barcode = clustScanpy[,1],
  umap_X = umapScanpy$X0,
  umap_Y = umapScanpy$X1,
  final_clusters = paste0("cluster_",clustScanpy$louvain)
)
rownames(meta) = as.character(meta$barcode)


barcodesPass = rownames(meta)[as.character(meta$final_clusters)%in%c("cluster_2","cluster_9","cluster_13")]

#Write data
write.table(barcodesPass, file = paste0(wd,"PhD_BPS53/data/revision1/18_endothelium_analysis/embryo_revision1_endothelium_erythroid_allantois.xgi"),
            append = FALSE,
            quote = FALSE, sep = "\t", eol = "\n",
            na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

cat("\n filtering matrix barcodes\n")

m = m[,as.character(barcodesPass)]

cat("\n filtering matrix by peaks in more than 0 cells\n")

m = m[rowSums(m)>0,]

cat("\n")
cat(dim(m))
cat("\n")

# Prepare data
print('Creating cisTopic object...')
count.matrix <- m

cisTopicObject <- createcisTopicObject(count.matrix = count.matrix, 
                                       project.name = "ATACseq",
                                       min.cells = 1,
                                       min.regions = 1,
                                       is.acc = 1,
                                       keepCountsMatrix=F)

cat("\n Running models\n")

cisTopicObject <- runModels(cisTopicObject, topic=c(50,60,70,80,90,100), 
                            seed=987, nCores=6, 
                            burnin = 120, 
                            iterations = 150, addModels=FALSE)

save(cisTopicObject,file=paste0(wd,"18_endothelium_analysis/data/cisTopicObject_50_100_afterReadsInPeaksQC_ery_endo_allantois.rda"))

cat("\n Selecting model\n")
pdf(paste0(wd,"18_endothelium_analysis/plots/cisTopic_modelSelect_50_100_afterReadsInPeaksQC_ery_endo_allantois.pdf"))
cisTopicObject <- selectModel(cisTopicObject)
dev.off()

cat("\n Log likelihood\n")

pdf(paste0(wd,"18_endothelium_analysis/plots/cisTopic_logLikelihoodByIter_50_100_afterReadsInPeaksQC_ery_endo_allantois.pdf"))

logLikelihoodByIter(cisTopicObject, select=c(50,60,70,80,90,100))
dev.off()

save(cisTopicObject,file=paste0(wd,"18_endothelium_analysis/data/cisTopicObject_50_100_afterReadsInPeaksQC_ery_endo_allantois.rda"))

cat("\n Convert to matrix\n")

cellassign <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
cellassign <- as(cellassign, "sparseMatrix")    

cat("\n Write matrix\n")

Matrix::writeMM(cellassign,
                file=paste0(wd,"18_endothelium_analysis/data/cisTopic_matrix_50_100_afterReadsInPeaksQC_ery_endo_allantois.mtx"))

