###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Output barcodes per cell type
##               "final_clusters" refers to clusters computed in
##               26_Final_visualisation_and_clustering.html
##               with louvain clustering.
###############################################################################

wd <- "/path/to/directory/"

meta = read.csv(paste0(wd,"PhD_BPS53/data/revision1/snATACseq_embryo_revision1_metadata_afterAllQC_passQC_UMAP_cluster.txt"),
                sep="\t")
rownames(meta)=meta$barcode
meta$final_clusters = paste0("cluster_",meta$final_clusters)

table(meta$final_clusters)

celltypes = unique(as.character(meta$final_clusters))
for (i in celltypes){
  bc = rownames(meta)[as.character(meta$final_clusters)==as.character(i)]
  i=  gsub("/","_",i)
  i=  gsub("\\.","_",i)
  i=  gsub(" ","_",i)

  write.table(bc,file=paste0(wd,"PhD_BPS53/data/revision1/17_BigWig_celltypes/",i,".txt"),
              sep='\n',quote=F,row.names=F,col.names=F)
}
