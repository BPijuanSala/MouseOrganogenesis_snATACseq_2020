###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  October 2019
##  DESCRIPTION: Compute TF enrichment scores in allantoic-haemato-endothelial landscape
###############################################################################

library(plyr)
suppressWarnings(library(cisTopic))
library(Matrix)

wd = "/path/to/directory/sample_pooled_preprocess_revision1/"
cat("\n loading\n")

load(file=paste0(wd,"18_endothelium_analysis/data/cisTopicObject_50_100_afterReadsInPeaksQC_ery_endo_allantois.rda"))
pred.matrix <- predictiveDistribution(cisTopicObject)


cat("\n loading peaks\n")

# Obtain signatures
path_to_signatures <- paste0(wd, '19_peaks_other_datasets/peaks_bed/')
Bulk_ATAC_signatures <- paste(path_to_signatures, list.files(path_to_signatures), sep='')
labels  <- gsub('_peaks_mm10.bed', '', list.files(path_to_signatures))

cat("\n getSignaturesRegions\n")

cisTopicObject <- getSignaturesRegions(cisTopicObject, 
                                       Bulk_ATAC_signatures, 
                                       labels=labels, 
                                       minOverlap = 0.4)

# To only keep unique peaks per signature
cisTopicObject@signatures <- llply(1:length(cisTopicObject@signatures), 
                                   function (i) cisTopicObject@signatures[[i]][-which(cisTopicObject@signatures[[i]] %in% unlist(as.vector(cisTopicObject@signatures[-i])))]) 
names(cisTopicObject@signatures) <- labels


cat("\n Compute cell rankings\n")

# Compute cell rankings 
library(AUCell)
aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)

cat("\n Check signature enrichment in cells\n")

# Check signature enrichment in cells 
cisTopicObject <- signatureCellEnrichment(cisTopicObject, aucellRankings, selected.signatures='all', aucMaxRank = 0.3*nrow(aucellRankings), plot=FALSE)

save(cisTopicObject,file=paste0(wd,"18_endothelium_analysis/data/cisTopicObject_50_100_afterReadsInPeaksQC_ery_endo_allantois_ChIPseq_peaks.rda"))


cellassign <- as(as.matrix(cisTopicObject@cell.data), "sparseMatrix")    


#Write results
Matrix::writeMM(cellassign,
                file=paste0(wd,"18_endothelium_analysis/data/cisTopic_matrix_50_100_afterReadsInPeaksQC_ery_endo_allantois_ChIPseq_enrichment_cell.mtx"))

write.table(rownames(cellassign),file=paste0(wd,"18_endothelium_analysis/data/cisTopic_matrix_50_100_afterReadsInPeaksQC_ery_endo_allantois_ChIPseq_enrichment_cell_cellNames.txt"),
            sep="\n",quote=F,row.names=F,col.names=F)

write.table(colnames(cellassign),file=paste0(wd,"18_endothelium_analysis/data/cisTopic_matrix_50_100_afterReadsInPeaksQC_ery_endo_allantois_ChIPseq_enrichment_cell_TF_colnames.txt"),
            sep="\n",quote=F,row.names=F,col.names=F)

pdf(paste0(wd,"18_endothelium_analysis/plots/cisTopic_modelSelect_50_100_afterReadsInPeaksQC_ery_endo_allantois_ChIPseq_heatmap.pdf"))

signaturesHeatmap(cisTopicObject)
dev.off()


####==========================
### LOCALLY
####==========================

all_colours = c(
  "Allantois" = "#532C8A",
  "Cardiomyocytes" =  "#B51D8D",  
  "Mixed mesoderm" =  "#ab80b7",#[26] E8.25 mixed mesoderm   
  "Endothelium" =  "#ff891c",                             
  "Erythroid" =  "#EF4E22",
  "ExE endoderm" = "#7F6874",                              
  "ExE mesoderm" =  "#8870ad",
  "Forebrain" = "#647a4f",
  "Mid/Hindbrain"="#8ca376",
  "Gut" =  "#EF5A9D",                                
  "Neural crest"= "#C3C388",
  "NMP" =  "#8EC792",                                        
  "Notochord" =  "#0F4A9C",                             
  "Paraxial mesoderm" =  "#8DB5CE",
  "Mesenchyme" = "#cc7818",
  "Somitic mesoderm" =  "#005579",                                   
  "Spinal cord" =  "#CDE088",                              
  "Surface ectoderm" = "#f7f79e",      
  "Pharyngeal mesoderm" =  "#C9EBFB"               
  
  
  
)

clustAnn_ATAC = c(
  
  "0"="NMP" ,
  "1"="Mixed mesoderm" , 
  "2"="Erythroid",
  
  '3'="Forebrain" , 
  "4"="Somitic mesoderm",
  "5"="Surface ectoderm",
  
  "6"="Paraxial mesoderm" ,
  
  "7"="Gut" ,  
  "8"="Spinal cord" ,          
  "9"="Endothelium"  ,
  "10"="Pharyngeal mesoderm" ,                             
  
  "11"="Mid/Hindbrain" ,
  "12"="Mesenchyme" ,
  "13"="Allantois",
  
  "14"="Cardiomyocytes" ,                                 
  "15"="Neural crest",
  
  "16"="ExE endoderm",
  "17"="Notochord"
)


wd <- "/path/to/directory/"

# Read files
m = Matrix::readMM(file=paste0(wd,"PhD_BPS53/data/revision1/18_endothelium_analysis/cisTopic_matrix_50_100_afterReadsInPeaksQC_ery_endo_allantois_ChIPseq_enrichment_cell.mtx"))

rownames(m) = read.table(file=paste0(wd,"PhD_BPS53/data/revision1/18_endothelium_analysis/cisTopic_matrix_50_100_afterReadsInPeaksQC_ery_endo_allantois_ChIPseq_enrichment_cell_cellNames.txt"),
            sep="\n",header=F)[,1]
colnames(m) = read.table(file=paste0(wd,"PhD_BPS53/data/revision1/18_endothelium_analysis/cisTopic_matrix_50_100_afterReadsInPeaksQC_ery_endo_allantois_ChIPseq_enrichment_cell_TF_colnames.txt"),
                         sep="\n",header=F)[,1]

#Read metadata
meta = read.table(paste0(wd,"PhD_BPS53/data/revision1/snATACseq_embryo_revision1_metadata_afterAllQC_passQC_UMAP_cluster.txt"),
                  sep="\t",header=T)
meta$final_clusters = paste0("cluster_",meta$final_clusters)
rownames(meta) = meta$barcode


#Plot ChIP-seq enrichment scores
bluePal <- c("#BFBFBF","#6495ED","#000000")

library(ggplot2)
colors = list()
for (i in 3:ncol(m)){
  chip=colnames(m)[i]
  df = data.frame(x = meta[rownames(m),"umap_X"],
                  y=meta[rownames(m),"umap_Y"],
                  value=m[,i])
  
  
  df = df[order(df$value,decreasing=F),]
  tiff(paste0(wd,"PhD_BPS53/plots/revision1/18_endothelium_analysis/UMAP_revision01_endothelium_ery_allantois_ChIPseq_",chip,"_absolute_withLegend.tiff"),
       width=7,height=7,res=300,units="in")
  
  p = ggplot(df, aes(x=x, y=y))
  plot = p + geom_point(aes(col = value))+
    theme_void() +
    scale_colour_gradientn(colours = bluePal
                           ,limits=c(0,1)) +
    #) +
    
    #theme(legend.title = element_blank(),legend.position = "none") + 
    ggtitle("")
  print(plot)

  dev.off()
  d = ggplot_build(plot)
  colors[[chip]] = (d$data[[1]][order(df$value,decreasing=F),"colour"])
  
}


# Compute quantiles of colours to make gradient of colours in PAGA 
# in python so that colours range from 0 to 1 in all plots

colors[["TAL1_HE"]][1]
colors[["TAL1_HE"]][length(colors[["TAL1_HE"]])/4]
colors[["TAL1_HE"]][length(colors[["TAL1_HE"]])/2]
colors[["TAL1_HE"]][length(colors[["TAL1_HE"]])*3/4]
colors[["TAL1_HE"]][length(colors[["TAL1_HE"]])]


colors[["TAL1_HP"]][1]
colors[["TAL1_HP"]][length(colors[["TAL1_HP"]])/4]
colors[["TAL1_HP"]][length(colors[["TAL1_HP"]])/2]
colors[["TAL1_HP"]][length(colors[["TAL1_HP"]])*3/4]
colors[["TAL1_HP"]][length(colors[["TAL1_HP"]])]


colors[["TAL1_Haemangioblast"]][1]
colors[["TAL1_Haemangioblast"]][length(colors[["TAL1_Haemangioblast"]])/4]
colors[["TAL1_Haemangioblast"]][length(colors[["TAL1_Haemangioblast"]])/2]
colors[["TAL1_Haemangioblast"]][length(colors[["TAL1_Haemangioblast"]])*3/4]
colors[["TAL1_Haemangioblast"]][length(colors[["TAL1_Haemangioblast"]])]

