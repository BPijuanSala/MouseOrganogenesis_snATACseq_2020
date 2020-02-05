###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to perform Quality Control based on reads in peaks
###############################################################################

wd <- "/path/to/directory/"

library(moduleColor)
library(ggplot2)
bluePal <- c("#BFBFBF","#6495ED","#000000")

#Read metadata
m = read.csv(file=paste0(wd,"data/revision1/07_doublet_removal/snATACseq_embryo_revision01_doublets_cisTopic_50_100_metadata.csv"),
             header=T)
rownames(m)<-m$index
m$doublet = rep(FALSE,nrow(m))

#Plot histogram with doublet threshold
m$doublet[m$doublet_scores > 0.4] = TRUE

hist((m$doublet_scores),breaks=50,col="black",main="",xlab="Doublet scores")
abline(v=log(0.4),col="red",lwd=2)



pdf(paste0(wd,"plots/revision1/QC_doublets.pdf"),width=5,height=4)

hist((m$doublet_scores),breaks=50,col="black",main="",xlab="Doublet scores")
abline(v=0.4,col="red",lwd=2)
dev.off()


#########Read data from first QC (this will read all barcodes)

constitutive_promoters <- read.table(paste0(wd,"data/mm10_consecutive_promoters.bed"))
num_of_reads = read.table(paste0(wd,"data/revision1/embryo_revision1.reads_per_cell"))
promoter_coverage = read.table(paste0(wd,"data/revision1/embryo_revision1.promoter_cov"))

qc = num_of_reads;
colnames(qc) <- c("barcode", "num_of_reads")
qc$promoter_coverage = 0;
#qc$read_in_peak = 0;
qc$read_in_promoter[match(promoter_coverage$V1, qc$barcode)] = promoter_coverage$V2
qc$ratio_promoters = qc$read_in_promoter/qc$num_of_reads

qc$promoter_coverage[match(promoter_coverage$V1, qc$barcode)] = promoter_coverage$V2/nrow(constitutive_promoters)

rownames(qc) = qc$barcode

m$barcode = m$index
m = m[,(colnames(m) %in% c("index","num_of_reads","promoter_coverage")) == F]
##########

#Merge metadata after first QC with data from first QC (this latter contains all barcodes)
meta_all = merge(m, qc, by="barcode",all=T)
############

######BARCODES
#Add information on the sorting gates during FACS

library(xlsx)
library(GenomicRanges)
bc1 <- as.character(read.xlsx(paste0(wd,"Experiment_data/Barcodes/Barcodes_plates.xlsx"),sheetName = 1)[,2])
bc1rev <- c()
for (i in bc1){ bc1rev <- c(bc1rev,reverse(chartr("ATGC","TACG",i)))}

bc2 <- as.character(read.xlsx(paste0(wd,"Experiment_data/Barcodes/Barcodes_plates.xlsx"),sheetName = "I2_E85_embryo_all")[,2])
bc2rev <- c()
for (i in bc2){ bc2rev <- c(bc2rev,reverse(chartr("ATGC","TACG",i)))}

bcSmall2 <-  as.character(read.xlsx(paste0(wd,"Experiment_data/Barcodes/Barcodes_plates.xlsx"),sheetName = "I2_E85_embryo_smallnuclei_2n")[,2])
bcSmall2Rev <- c()
for (i in bcSmall2){bcSmall2Rev <- c(bcSmall2Rev,reverse(chartr("ATGC","TACG",i)))}

bcBig2 <-  as.character(read.xlsx(paste0(wd,"Experiment_data/Barcodes/Barcodes_plates.xlsx"),sheetName = "E85_embryo_bignuclei_4n")[,2])
bcBig2Rev <- c()
for (i in bcBig2){bcBig2Rev <- c(bcBig2Rev,reverse(chartr("ATGC","TACG",i)))}

qc$nuclei_type = rep("Unassigned",nrow(qc))

meta_all[substr(as.character(meta_all$barcode),11,22) %in% bcSmall2Rev,"nuclei_type"] = "2n"
meta_all[substr(as.character(meta_all$barcode),11,22) %in% bcBig2Rev,"nuclei_type"] = "4n"
meta_all[substr(as.character(meta_all$barcode),11,22) %in% bc2rev,"nuclei_type"] = "all"

table(meta_all$nuclei_type,meta_all$doublet,useNA="always")


###########
####Add information of number of reads in peaks

############

read_in_peak = read.table(paste0(wd,"data/revision1/12_barcodeStats_celltypePeaks/embryo_celltypeSpecificPeaks.reads_in_peak"))

meta_all$read_in_peak[match(read_in_peak$V1, meta_all$barcode)] = read_in_peak$V2
meta_all$ratio_peaks = as.numeric(meta_all$read_in_peak)/meta_all$num_of_reads


#Filter data that passed the first QC
meta = meta_all[(meta_all$doublet %in% c(TRUE,NA)) == F,]
dim(meta)

rownames(meta) =meta$barcode

#Add UMAP computed after QC and clusters
umapScanpy = read.csv(paste0(wd,"data/revision1/07_doublet_removal/snATACseq_embryo_revision1_firstQC_afterDoubletRemoval_UMAP.csv"))
clustScanpy = read.csv(paste0(wd,"data/revision1/07_doublet_removal/snATACseq_embryo_revision1_firstQC_afterDoubletRemoval_Louvain.csv"))
rownames(umapScanpy)<-clustScanpy[,1]
meta$umap_X = meta$umap_Y = NA
meta[rownames(umapScanpy),"umap_X"] = umapScanpy[,2]
meta[rownames(umapScanpy),"umap_Y"] = umapScanpy[,3]

#Plot histograms with thresholds
hist(log(meta$num_of_reads+1),breaks=50,col="black",main="",xlab="log(number of reads+1)")
abline(v=log(2000+1),col="red")

hist(log(meta$promoter_coverage),breaks=100,col="black",main="",xlab="Reads in constitutive promoters / Total constitutive promoters")
abline(v=log(0.03),col="red")

par(mfrow=c(1,1))
hist(log(meta$ratio_peaks),col="black",breaks=100,main="",xlab="Ratio reads in peaks")
abline(v=log(0.24),col="red",lwd=3)


pdf(paste0(wd,"plots/revision1/QC_ratioPeaks.pdf"),width=5,height=4)
hist(log(meta$ratio_peaks),breaks=50,col="black",main="",xlab="log(ratio of peaks)")
abline(v=log(0.24),col="red",lwd=2)
dev.off()


#Set ratio peaks threshold to 0.24
thres=0.24
metaNew1 = meta[meta$ratio_peaks<thres,]

#Plot nuclei before filtering by second QC
plot(as.numeric(as.character(meta$umap_X)),
     as.numeric(as.character(meta$umap_Y)),pch=20,col="black",cex=0.7)
points(as.numeric(as.character(metaNew1$umap_X)),
     as.numeric(as.character(metaNew1$umap_Y)),pch=20,col="red",cex=0.7)

#Plot nuclei after filtering by second QC

metaNew = meta[meta$ratio_peaks>thres,]
plot(as.numeric(as.character(metaNew$umap_X)),
     as.numeric(as.character(metaNew$umap_Y)),pch=20,col="black",cex=0.7)

#Plot ratio_peaks on nuclei after filtering by second QC

p = ggplot(metaNew, aes(x=umap_X, y=umap_Y))
plot = p + geom_point(aes(col = ratio_peaks))+
  theme_void() +
  scale_colour_gradientn(colours = bluePal,trans="log2") 
print(plot)

#Plot doublet_scores on nuclei after filtering by second QC

p = ggplot(metaNew, aes(x=umap_X, y=umap_Y))
plot = p + geom_point(aes(col = doublet_scores))+
  theme_void() +
  scale_colour_gradientn(colours = bluePal,trans="log2") 
print(plot)



#Write data
#Output nuclei barcodes that passed QC
write.table(rownames(meta)[meta$ratio_peaks>0.24], file = paste0(wd,"data/revision1/12_barcodeStats_celltypePeaks/embryo_revision1_readsPeaks24.xgi"),
            append = FALSE,
            quote = FALSE, sep = "\t", eol = "\n",
            na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

#write metadata of cells that passed QC

write.table(metaNew,paste0(wd,"data/snATACseq_embryoTogether_allPeaks_metadata_round4_afterDoubletRemoval_cellTypePeakCall_filter.txt"),
            quote=F,sep="\t",col.names=T,row.names=F)


meta_all$pass_all_QC = rep(FALSE,nrow(meta_all))
meta_all[match(rownames(meta)[meta$ratio_peaks>0.24],meta_all$barcode),"pass_all_QC"] = TRUE


#Write metadata of all nuclear barcodes
write.table(meta_all,paste0(wd,"data/revision1/snATACseq_embryo_revision1_metadata_afterAllQC.txt"),
            quote=F,sep="\t",col.names=T,row.names=F)



