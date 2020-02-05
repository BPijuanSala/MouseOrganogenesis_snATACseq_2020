###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to assess quality of nuclei and discard those
##               of low quality.
###############################################################################


###Script adapted from from https://github.com/r3fang/snATAC

##################################################################
# 1) reads per barcode > 2000
# 2) constitutive_promoters promoter coverage > 3%
##################################################################

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

#Plotting
pdf(paste0(wd,"plots/revision1/QC_numReads.pdf"),width=5,height=4)
hist(log(qc$num_of_reads+1),breaks=50,col="black",main="",
     xlab="log(number of reads+1)")
abline(v=log(2000+1),col="red")
dev.off()

pdf(paste0(wd,"plots/revision1/QC_promCov.pdf"),width=5,height=4)

hist(log(qc$promoter_cov),breaks=100,col="black",main="",xlab="Reads in constitutive promoters / Total constitutive promoters")
abline(v=log(0.03),col="red")
dev.off()

tiff(paste0(wd,"plots/revision1/QC_promCov_numReads.tiff"),width=5,height=4,
     res=200,units="in")

plot(log(qc$num_of_reads),log(qc$promoter_cov),pch=20,cex=0.2,#ylab="log(Promoter coverage)",
     #xlab="log(Number of reads)",
     xlim=c(0,13),ylim=c(-9,0),xlab="",ylab="",
     xaxt="n",yaxt="n")
#rect(xleft=log(2000),xright=13,
#          ybottom=log(0.03),ytop=0.11,border="red")
dev.off()


pdf(paste0(wd,"plots/revision1/QC_promCov_numReads_frame.pdf"),width=5,height=4
     )

plot(log(qc$num_of_reads)[1],log(qc$promoter_cov)[1],pch=20,cex=0.2,
     #ylab="log(Promoter coverage)",
     #xlab="log(Number of reads)",
     xlim=c(0,13),ylim=c(-9,0),xlab="",ylab=""
     #xaxt="n",yaxt="n"
     )
rect(xleft=log(2000),xright=13,
          ybottom=log(0.03),ytop=0.11,border="red")
dev.off()

rownames(qc) =qc$barcode



##################
#Set thresholds
idx <- which(qc$promoter_cov > 0.03 & qc$num_of_reads > 2000)


length(idx)

#Select cells that pass QC.
qc_sel <- qc[idx,]


#Write data
write.table(qc_sel[,1], file = paste0(wd,"data/revision1/embryo_revision1.xgi"),
            append = FALSE,
            quote = FALSE, sep = "\t", eol = "\n",
            na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
