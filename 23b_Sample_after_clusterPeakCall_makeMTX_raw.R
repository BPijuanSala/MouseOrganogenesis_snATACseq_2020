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
mat1 <- Matrix::readMM(paste0(wd,"11_matrix_afterClusterQC/embryo_revision1_allPeaks_afterclusterPeak_passedQC_subset00_01_raw.mtx"))
cat("done first matrix \n")
cat(dim(mat1))

mat2 <- Matrix::readMM(paste0(wd,"11_matrix_afterClusterQC/embryo_revision1_allPeaks_afterclusterPeak_passedQC_subset02_03_raw.mtx"))

mat3 = Matrix::rbind2(mat1,mat2)

Matrix::writeMM(mat3,file=paste0(wd,"11_matrix_afterClusterQC/embryo_revision1_allPeaks_afterClusterPeak_raw.mtx"))
cat("DONE third script \n")
