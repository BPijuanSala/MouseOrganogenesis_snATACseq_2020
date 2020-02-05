###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to run cisTopic on the final amtrix
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

barcodesPass = as.vector(read.csv(paste0(wd,"12_barcodeStats_celltypePeaks/embryo_revision1_readsPeaks24.xgi"),
                                  sep="\n",header=F)[,1])  

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

save(cisTopicObject,file=paste0(wd,"13_cisTopic/data/cisTopicObject_50_100_afterReadsInPeaksQC_24.rda"))

cat("\n Selecting model\n")
pdf(paste0(wd,"13_cisTopic/plots/cisTopic_modelSelect_50_100_afterReadsInPeaksQC_24.pdf"))
cisTopicObject <- selectModel(cisTopicObject)
dev.off()
cat("\n Log likelihood\n")

pdf(paste0(wd,"13_cisTopic/plots/cisTopic_logLikelihoodByIter_50_100_afterReadsInPeaksQC_24.pdf"))

logLikelihoodByIter(cisTopicObject, select=c(50,60,70,80,90,100))
dev.off()

save(cisTopicObject,file=paste0(wd,"13_cisTopic/data/cisTopicObject_50_100_afterReadsInPeaksQC_24.rda"))

cat("\n Convert to matrix\n")

cellassign <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
cellassign <- as(cellassign, "sparseMatrix")    

cat("\n Write matrix\n")

Matrix::writeMM(cellassign,
                file=paste0(wd,"13_cisTopic/data/cisTopic_matrix_50_100_afterReadsInPeaksQC_24.mtx"))

