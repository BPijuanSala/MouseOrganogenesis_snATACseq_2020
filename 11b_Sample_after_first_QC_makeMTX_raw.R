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
mat1 <- Matrix::readMM(paste0(wd,"06_matrix/embryo_revision1_allPeaks_passedQC_subset00_01_raw.mtx"))
cat("\n\n done first matrix \n")
cat(dim(mat1))

mat2 <- Matrix::readMM(paste0(wd,"06_matrix/embryo_revision1_allPeaks_passedQC_subset02_03_raw.mtx"))
cat("\n\n done second matrix \n")

cat(dim(mat2))

mat3 = Matrix::rbind2(mat1,mat2)
cat(dim(mat3))

Matrix::writeMM(mat3,file=paste0(wd,"06_matrix/embryo_revision1_allPeaks_passedQC_raw.mtx"))
cat("DONE third script \n")
