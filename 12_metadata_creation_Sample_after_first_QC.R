###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to create metadata.
###############################################################################


wd <- "/path/to/directory/"


constitutive_promoters <- read.table(paste0(wd,"data/mm10_consecutive_promoters.bed"))#This is a supplementary file in the publication.
num_of_reads = read.table(paste0(wd,"data/revision1/embryo_revision1.reads_per_cell"))
promoter_cov = read.table(paste0(wd,"data/revision1/embryo_revision1.promoter_cov"))

qc = num_of_reads; 
colnames(qc) <- c("barcode", "num_of_reads")
qc$promoter_cov = 0; 

qc$read_in_promoter[match(promoter_cov$V1, qc$barcode)] = promoter_cov$V2
qc$ratio_promoters = qc$read_in_promoter/qc$num_of_reads

qc$promoter_cov[match(promoter_cov$V1, qc$barcode)] = promoter_cov$V2/nrow(constitutive_promoters)



meta = qc


idx <- which(meta$promoter_cov > 0.03 & meta$num_of_reads > 2000)

length(idx)
qc_sel <- qc[idx,]


write.table(qc_sel,paste0(wd,"data/revision1/metadata_revision1_firstQC.txt"),
            quote=F,sep="\t",col.names=T,row.names=F)




