###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Convert matrix to mtx
###############################################################################

library(Matrix)
library(matrixStats)
library(irlba)



wd <- "/path/to/directory/sample_pooled_preprocess_revision1/"
cat("reading \n")
mat1 <- t(as.matrix(read.table(paste0(wd,"11_matrix_afterClusterQC/embryo_revision1_allPeaks_afterclusterPeak_passedQC_subset00.mat.bin"))))
cat("done first matrix \n")
cat(dim(mat1))

matS1 <- as(mat1, "sparseMatrix")
rm(mat1)

mat2 <- t(as.matrix(read.table(paste0(wd,"11_matrix_afterClusterQC/embryo_revision1_allPeaks_afterclusterPeak_passedQC_subset01.mat.bin"))))
cat("done second matrix \n")
cat(dim(mat2))

matS2 <- as(mat2, "sparseMatrix")
rm(mat2)

matS3 = Matrix::rbind2(matS1,matS2)

Matrix::writeMM(matS3,file=paste0(wd,"11_matrix_afterClusterQC/embryo_revision1_allPeaks_afterclusterPeak_passedQC_subset00_01.mat.bin.mtx"))

cat("DONE first script \n")
