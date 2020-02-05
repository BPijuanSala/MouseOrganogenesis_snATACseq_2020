###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  December 2019
##  DESCRIPTION: Run cisTopic on the first binary matrix
###############################################################################


############################################
# Install packages if you don't have them
############################################

#devtools::install_github("aertslab/AUCell")
#devtools::install_github("aertslab/RcisTarget")
#source("https://bioconductor.org/biocLite.R")
#biocLite('GenomicRanges')

#devtools::install_github("aertslab/cisTopic")
#source("https://bioconductor.org/biocLite.R")
#biocLite(c('Rsubread', 'umap', 'Rtsne', 'ComplexHeatmap', 'fastcluster', 'data.table', 'rGREAT', 'ChIPseeker', 'TxDb.Hsapiens.UCSC.hg19.knownGene', 'org.Hs.eg.db', 'densityClust'))


#Load libraries
suppressWarnings(library(cisTopic))
library(Matrix)


wd = "/path/to/directory/sample_pooled_preprocess_revision1/"


cat("\n reading matrix\n")
m <- Matrix::readMM(paste0(wd, '06_matrix/embryo_revision1_allPeaks_passedQC.mat.bin.mtx'))*1


cat(dim(m))

cat("\n reading barcodes\n")

barcodes = as.vector(read.csv(paste0(wd,"06_matrix/embryo_revision1_allPeaks_passedQC_barcodeNames.xgi"),
                    sep="\n",header=F)[,1])   

cat(dim(barcodes))

colnames(m) <- barcodes

cat("\n reading bed file\n")

regions_frame = read.table(paste0(wd,"05_MACS2_afterQC/embryo_revision1_nuclearGenes_passedQC_summits_extended500bp_overlapped_noblacklist_TSS_merged.bed"),
                                  header=F, sep="\t")
colnames(regions_frame) <- c('Chr', 'Start', 'End')

rownames(m) <- paste0(regions_frame[,1], ':', regions_frame[,2], '-', regions_frame[,3])

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

cisTopicObject <- runModels(cisTopicObject, topic=c(20,30, 40, 50,60,70), 
                            seed=987, nCores=6, 
                            burnin = 120, 
                            iterations = 150, addModels=FALSE)

save(cisTopicObject,file=paste0(wd,"07_doublet_removal/data/cisTopicObject.rda"))

cat("\n Selecting model\n")
pdf(paste0(wd,"07_doublet_removal/plots/cisTopic_modelSelect.pdf"))
cisTopicObject <- selectModel(cisTopicObject)
dev.off()
cat("\n Log likelihood\n")

pdf(paste0(wd,"07_doublet_removal/plots/cisTopic_logLikelihoodByIter.pdf"))

logLikelihoodByIter(cisTopicObject, select=c(20,30, 40, 50,60,70))
dev.off()

save(cisTopicObject,file=paste0(wd,"07_doublet_removal/data/cisTopicObject.rda"))

cat("\n Convert to matrix\n")

cellassign <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
cellassign <- as(cellassign, "sparseMatrix")    

cat("\n Write matrix\n")

Matrix::writeMM(cellassign,
                file=paste0(wd,"07_doublet_removal/data/cisTopic_matrix.mtx"))


