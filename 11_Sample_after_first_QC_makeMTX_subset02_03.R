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
mat1 <- t(as.matrix(read.table(paste0(wd,"06_matrix/embryo_revision1_allPeaks_passedQC_subset02.mat.bin"))))
cat("done first matrix \n")
cat(dim(mat1))

matS1 <- as(mat1, "sparseMatrix")
rm(mat1)

mat2 <- t(as.matrix(read.table(paste0(wd,"06_matrix/embryo_revision1_allPeaks_passedQC_subset03.mat.bin"))))
cat("done second matrix \n")
cat(dim(mat2))

matS2 <- as(mat2, "sparseMatrix")
rm(mat2)

matS3 = Matrix::rbind2(matS1,matS2)

Matrix::writeMM(matS3,file=paste0(wd,"06_matrix/embryo_revision1_allPeaks_passedQC_subset02_03.mat.bin.mtx"))
cat("DONE second script \n")