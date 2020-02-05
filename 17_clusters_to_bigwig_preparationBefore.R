###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Output barcodes per cell type
##               "ann2" refers to clusters computed in
##               16_Visualisation_clustering_after_doublet_removal.html
##               with louvain clustering.
###############################################################################

wd <- "/path/to/directory/"

meta = read.csv(file=paste0(wd,"data/revision1/07_doublet_removal/snATACseq_embryo_revision01_doublets_cisTopic_50_100_metadata_afterDoubletRemoval.csv"))
rownames(meta)<-meta$index

meta$louvain_afterDoubletRemoval = paste0("cluster_",meta$louvain_afterDoubletRemoval)

celltypes = unique(as.character(meta$louvain_afterDoubletRemoval))
for (i in celltypes){
  bc = rownames(meta)[as.character(meta$louvain_afterDoubletRemoval)==as.character(i)]
  i=  gsub("/","_",i)
  i=  gsub("\\.","_",i)
  i=  gsub(" ","_",i)

  write.table(bc,file=paste0(wd,"data/revision1/cluster_barcode_files/",i,".txt"),
              sep='\n',quote=F,row.names=F,col.names=F)
}

