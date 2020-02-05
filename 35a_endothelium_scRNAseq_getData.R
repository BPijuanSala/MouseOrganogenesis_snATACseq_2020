
###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  October 2019
##  DESCRIPTION: Get allantois, mixed mesoderm, endothelium, haemato-endothelium
###############################################################################


library(ggplot2)
library(Matrix)
library(knitr)
library(scran)
library(igraph)
library(Rtsne)
wd <- "/path/to/directory/scRNAseq_PijuanSala_2019/"
dir = "/path/to/directory/sample_pooled_preprocess_revision1/"
#Get data
meta = read.csv(paste0(wd,"data/PijuanSalaEtAl_SupplementaryTable3_revision.txt"), 
                  header = T, stringsAsFactors = F, sep = "\t")
meta0 = read.table(paste0(wd,"data/rawData/meta.tab"),header=T,stringsAsFactors=F,sep="\t")
counts = readMM(paste0(wd,"data/rawData/raw_counts.mtx"))

# dim(counts)
#[1]  27998 127581
counts = as(counts, "dgCMatrix")


sfs = as.vector(read.table(paste0(wd,"data/rawData/sizefactors.tab"))[,1])
names(sfs) <- meta0$cell
genes = read.table(paste0(wd,"data/rawData/genes.tsv"), stringsAsFactors = FALSE)
colnames(counts) = meta0$cell
rownames(counts) = genes[,1]


rownames(meta) <- as.character(meta$cell)

##Select cells
cellsSel <- rownames(meta)[as.character(meta$celltype)%in%c("Mixed mesoderm","Endothelium","Haematoendothelial progenitors","Allantois")]
countsSel <- counts[,cellsSel]
sfsSel <- sfs[cellsSel]


#Normalise subset
library(scran)
#library(scmap)
library(SingleCellExperiment)

sce = SingleCellExperiment(assays = list("counts" = countsSel))
sizeFactors(sce) = sfsSel
sce = scater::normalize(sce)



save(sce,file=paste0(dir,"21_endothelium_scRNAseq/data/sce_norm_endothelium.rda"))
#load(paste0(wd,"data/sce_norm.rda"))


###Save data

countsData <- as(round(SingleCellExperiment::logcounts(sce),2),"sparseMatrix")

writeMM(countsData,file=paste0(dir,"21_endothelium_scRNAseq/data/countsData_endothelium_logged.mtx"))
write.table(rownames(countsData),file=paste0(dir,"21_endothelium_scRNAseq/data/genes_endothelium.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\n")

write.table(colnames(countsData),file=paste0(dir,"21_endothelium_scRNAseq/data/cells_endothelium.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\n")