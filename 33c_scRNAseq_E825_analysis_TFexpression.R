###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to calculate the expression levels of the top 20 TFs from
##                Ibarra-Soria et al., Nat. Cell Biol., 2018 and 
##                Pijuan-Sala et al., Nature, 2019
###############################################################################

wd <- "/path/to/file/"
####====================
# Load data and normalise
####====================

meta <- read.csv(paste0(wd,"PhD_BPS32/release6/paper_resources/important_data/PijuanSalaEtAl_SupplementaryTable3_revision.txt"), 
                  sep = "\t")    
rownames(meta) <- meta$cell


counts <- Matrix::readMM(file=paste0(wd,"PhD_BPS53/data/05_counts_scRNAseq/counts_raw_E825.mtx"))

colnames(counts) <- as.character(read.table(paste0(wd,"PhD_BPS53/data/05_counts_scRNAseq/cells_E825.tab"),header=T)[,1])
rownames(counts) <-  as.character(read.table(paste0(wd,"PhD_BPS53/data/05_counts_scRNAseq/genes_E825.tab"),header=TRUE)[,1])
genesTab <-  read.table(paste0(wd,"PhD_BPS32/release6/data/genes.tsv"),header=F)
colnames(genesTab) = c("Gene.ID","Associated.Gene.Name")

meta <- meta[colnames(counts),]
genesRNAseq<-  as.character(read.table(paste0(wd,"PhD_BPS53/data/05_counts_scRNAseq/HVG_scRNAseq_E825_scran.tab"),header=TRUE)[,1])

library(SingleCellExperiment)

sce <- SingleCellExperiment::SingleCellExperiment(list(counts=counts))

lib.sizes = Matrix::colSums(counts(sce))
sce2 = sce[scater::calcAverage(sce)>0.1,]
dim(sce)
dim(sce2)
clusts = as.numeric(scran::quickCluster(sce2, method = "igraph", min.size = 100))

min.clust = min(table(clusts))/2

new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce2 = scran::computeSumFactors(sce2, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

sizeFactors(sce) = sizeFactors(sce2)

sce = normalize(sce)

require(biomaRt)
min.mean=10^-4
trend = scran::trendVar(sce, use.spikes = FALSE, loess.args = list(span = 0.05))
decomp = scran::decomposeVar(sce, fit = trend)
decomp = decomp[decomp$mean > min.mean,]


#exclude sex genes
xist = "ENSMUSG00000086503"

ychr = read.table(paste0(wd,"bin_data/mmusculus_Y_genes.txt"), 
                  stringsAsFactors = FALSE)[,1]
decomp = decomp[!rownames(decomp) %in% c(xist, ychr),]

decomp$FDR = p.adjust(decomp$p.value, method = "fdr")
HVG <- rownames(decomp)[decomp$p.value < 0.05]

length(HVG)
countsNorm = SingleCellExperiment::logcounts(sce)
#save(countsNorm,file=paste0(wd,"PhD_BPS53/data/countsNorm_E825_RNAseq10X.rda"))

###########==================
#load(file=paste0(wd,"PhD_BPS53/data/countsNorm_E825_RNAseq10X.rda"))

#Load the top TFs and explore their expression profiles
load(file=paste0(wd,"PhD_BPS53/data/revision1/16_chromVAR/snATACseq_enriched_TFs_top15_wilcoxon_allPeaks.rda"))
TFenrichedID = anSeq::getGeneID(TFenriched,Table = genesTab)[[2]]
TFenrichedName = anSeq::getGeneName(TFenrichedID,Table = genesTab)[[2]]
TFenriched[TFenriched%in%TFenrichedName ==F]
TFenrichedName[TFenrichedID%in%rownames(countsNorm)]
TFenrichedID[!TFenrichedID%in%rownames(countsNorm)]


countsSel=countsNorm[TFenrichedID,]
c = meta[colnames(countsSel),"celltype"]
heatCountsMean0 <- t(apply(countsSel,1,function(x){
  
  y <- aggregate(x,by=list(c),FUN=mean)
  y2 <- y$x[order(as.character(y$Group.1))]
  names(y2) <- as.character(y$Group.1)[order(as.character(y$Group.1))]
  
  
  return(y2)
}))

rownames(heatCountsMean0) = TFenrichedName
geneOrder = den$gene.clust$labels[den$gene.clust$order]
notThere = TFenriched[TFenriched%in%TFenrichedName ==F]

geneOrder = geneOrder[(geneOrder %in% notThere) ==F]
cellTypeOrder =den$cell.clust$labels[den$cell.clust$order]
cellTypeOrder = c("Notochord","NMP","Somitic mesoderm","Intermediate mesoderm","ExE mesoderm",
                  "Allantois","Endothelium","Erythroid1","Erythroid2","Erythroid3",
                  "ExE endoderm","Gut","Mesenchyme","Cardiomyocytes",
                  "Pharyngeal mesoderm","Paraxial mesoderm","Spinal cord",
                  "Forebrain/Midbrain/Hindbrain","Neural crest","Surface ectoderm"
                 )
heatCountsStd <- t(apply(heatCountsMean0[geneOrder,cellTypeOrder], 1, function(x) x/(max(x)+0.00001)))  # standarise
heatCountsStdRNA = heatCountsStd
heatCountsMeanRNA =heatCountsMean0[geneOrder,cellTypeOrder]

#Save the expression profiles
save(heatCountsStdRNA,heatCountsMeanRNA,file=paste0(wd,"PhD_BPS53/data/revision1/16_chromVAR/ExpressionTFs_RNA_celltypes_wilcoxon_20190310_peaksAll.rda"))

#PLEASE RETURN TO PREVIOUS SCRIPT!