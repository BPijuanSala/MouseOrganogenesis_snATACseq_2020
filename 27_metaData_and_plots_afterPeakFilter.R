###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to generate final metadata
###############################################################################


wd <- "/path/to/directory/"
library(moduleColor)
library(ggplot2)
library(anSeq)
bluePal <- c("#BFBFBF","#6495ED","#000000")

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


#Read metadata
meta_all = read.table(file=paste0(wd,"PhD_BPS53/data/revision1/snATACseq_embryo_revision1_metadata_afterAllQC.txt"),
             header=T,sep="\t")


#Filter metadata for cells that passed QC
m = meta_all[meta_all$pass_all_QC,]
rownames(m)<-as.character(m$barcode)
dim(m)

#Add final UMAP and clustering
umapScanpy<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/13_cisTopic/snATACseq_embryo_revision1_reads_in_peaks_above_24_UMAP.csv"))
clustScanpy<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/13_cisTopic/snATACseq_embryo_revision1_reads_in_peaks_above_24_Louvain.csv"))
rownames(umapScanpy)<-clustScanpy[,1]
m$final_clusters = m$umap_Y = m$umap_X = "None"
m[rownames(umapScanpy),"umap_X"] = umapScanpy$X0
m[rownames(umapScanpy),"umap_Y"] = umapScanpy$X1
m[rownames(umapScanpy),"final_clusters"] = clustScanpy$louvain


######ANNOTATION

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

#Add annotation to clustering
m$ann = clustAnn_ATAC[as.character(m$final_clusters)]


all_colours_sub = all_colours[clustAnn_ATAC[unique(as.character(m$final_clusters))]]

meta=m
meta$umap_X =as.numeric(as.character(meta$umap_X))
meta$umap_Y =as.numeric(as.character(meta$umap_Y))

meta$umap_X = -meta$umap_X



#Add metadata to metadata containing all barcodes
meta_all$ann = meta_all$final_clusters = meta_all$umap_Y = meta_all$umap_X = "None"
meta_all[match(as.character(meta$barcode),as.character(meta_all$barcode)),"ann"] = meta$ann
meta_all[match(as.character(meta$barcode),as.character(meta_all$barcode)),"final_clusters"] = meta$final_clusters
meta_all[match(as.character(meta$barcode),as.character(meta_all$barcode)),"umap_Y"] = meta$umap_Y
meta_all[match(as.character(meta$barcode),as.character(meta_all$barcode)),"umap_X"] = meta$umap_X


#Cleanup metadata

m = meta[,c('barcode','nuclei_type','num_of_reads',
            'promoter_coverage','read_in_promoter',
            'doublet_scores','read_in_peak','ratio_peaks',
            'umap_X','umap_Y','final_clusters','ann')]


meta_all2 = meta_all[,c('barcode','nuclei_type','num_of_reads',
            'promoter_coverage','read_in_promoter',
            'doublet_scores','doublet','read_in_peak','ratio_peaks', 'pass_all_QC',
            'umap_X','umap_Y','final_clusters','ann')]

#write.table(m,paste0(wd,"PhD_BPS53/data/revision1/Supplementary_Fig5_metadata_passQC.txt"),
 #      quote=F,sep="\t",col.names=T,row.names=F)
#write.table(m,paste0(wd,"PhD_BPS53/data/revision1/snATACseq_metadata_revision2_allNuclei.txt"),
 #      quote=F,sep="\t",col.names=T,row.names=F)


#Plot UMAP with all cell types
pdf(paste0(wd,"PhD_BPS53/plots/revision1/UMAP_revision01_celltypes.pdf"),
   width=8,height=7)
p = ggplot(meta, aes(x=umap_X, y=umap_Y))
p + geom_point(aes(col = factor(ann))) + 
  scale_color_manual(values = all_colours_sub[order(names(all_colours_sub))], 
                     labels = names(all_colours_sub)[order(names(all_colours_sub))], 
                     drop = FALSE, name = "")  +
  theme_void() +
  theme(legend.text=element_text(size=rel(1.1))
  ) + ggtitle("")
dev.off()


tiff(paste0(wd,"PhD_BPS53/plots/revision1/UMAP_revision01_clusters_celltypes.tiff"),
    width=7,height=7,res=600,units="in")
plot(meta$umap_X,meta$umap_Y,pch=20,cex=0.5,axes = F,ylab="",xlab="",
     col=all_colours[as.character(meta$ann)])

dev.off()


###TOPICS
#Plot topics in UMAP. e.g. topic 20

cistopic = Matrix::readMM(paste0(wd,"PhD_BPS53/data/revision1/13_cisTopic/cisTopic_matrix_50_100_afterReadsInPeaksQC_24.mtx"))
names_cist = as.character(read.csv(paste0(wd,"PhD_BPS53/data/revision1/12_barcodeStats_celltypePeaks/embryo_revision1_readsPeaks24.xgi"),
                      header=F)[,1])

colnames(cistopic) = names_cist
library(Matrix)

df = data.frame(x = meta[colnames(cistopic),"umap_X"],
                y=meta[colnames(cistopic),"umap_Y"],
                value=cistopic[20,])
df = df[order(df$value,decreasing=F),]

tiff(paste0(wd,"PhD_BPS53/plots/revision1/UMAP_revision01_topic20.tiff"),
     width=7,height=7,res=300,units="in")

p = ggplot(df, aes(x=x, y=y))
plot = p + geom_point(aes(col = value))+
  theme_void() +
  scale_colour_gradientn(colours = bluePal) +
  theme(legend.title = element_blank(),legend.position = "none") + 
  ggtitle("")
print(plot)
dev.off()

##==========================NUCLEI SIZE

#Plot data sorted in 2n and 4n gates
pdf(paste0(wd,"PhD_BPS53/plots/revision1/snATACseq_revision1_2n_4n_UMAP.pdf"),
    width=10,height=7)
par(mfrow=c(1,2))
plot(as.numeric(as.character(m$umap_X)),yaxt="n",xaxt="n",ylab="",xlab="",main="2n",
     as.numeric(as.character(m$umap_Y)),pch=20,col="gray",cex=0.7)

points(as.numeric(as.character(m$umap_X))[as.character(m$nuclei_type)=="2n"],
       as.numeric(as.character(m$umap_Y))[as.character(m$nuclei_type)=="2n"],
       col="red",
       pch=20,cex=0.7)

plot(as.numeric(as.character(m$umap_X)),yaxt="n",xaxt="n",ylab="",xlab="",main="4n",
     as.numeric(as.character(m$umap_Y)),pch=20,col="gray",cex=0.7)

points(as.numeric(as.character(m$umap_X))[as.character(m$nuclei_type)=="4n"],
       as.numeric(as.character(m$umap_Y))[as.character(m$nuclei_type)=="4n"],
       col="blue",
       pch=20,cex=0.7)
dev.off()

####NUCLEI SIZE
meta=m
metaSub = meta[as.character(meta$nuclei_type)%in%c("2n","4n"),]
metaSub$nuclei_type=as.character(metaSub$nuclei_type)
tab =table(metaSub$nuclei_type,metaSub$ann)
tab2 = apply(tab,2,function(x){x/sum(x)})


pdf(paste0(wd,"PhD_BPS53/plots/revision1/UMAP_all_together_clusters_afterDoubletRemoval_nuclei.pdf"),
    width=10,height=5)
par(mar=c(15,3,3,3))
barplot(tab2[c(2,1),order(tab2["4n",])],las=2)
dev.off()





############################
# Heatmap TSS
###########################
#Plot accessibility at known TSS used to annotate cell types

#Read total accessibility counts per region and per cell type
counts0 = read.table(paste0(wd,"PhD_BPS53/data/revision1/snATACseq_embryo_revision01_pooledNumber_OCR_perCelltype.csv"),
                     sep=",",header=T)
counts = counts0[,-1]
rownames(counts) = as.character(read.table(paste0(wd,"PhD_BPS53/data/revision1/snATACseq_embryo_revision01_pooledNumber_OCR_perCelltype_peakLabel.csv"),
                                           sep=",",header=T)[,2])
colnames(counts) = as.character(read.table(paste0(wd,"PhD_BPS53/data/revision1/snATACseq_embryo_revision01_pooledNumber_OCR_perCelltype_celltypeLabel.csv"),
                                           sep=",",header=T)[,2])

library(anSeq)#From https://github.com/BPijuanSala/anSeq

#Define list of genes to assess
genes = c(
  "Tfap2b","Sox10","Foxd3",#NC
  "Hbb-y","Hba-a2","Hba-a1","Gypa",#Erythroid
  "Cdh5","Pecam1","Kdr","Etv2", #Endothelium
  "Slc2a2","Fxyd2","Aqp3","Aldob",#ExE endoderm
  "Foxa2","T","Chrd","Noto","Shh",#,#Notochord
  "Foxc2","Col26a1","Meox1","Foxd2",#Paraxial
  "Tdo2",#Mesenchyme
  "Sdpr","Wisp1",#Mesenchyme
  "Ahnak", 
  "Nkx2-5","Tnnt2","Myl7",#Cardiomyocytes
  "Lor",#Intermediate meso
  "Hoxd1","Haglr",#ExE meso
  "Hoxa10","Tbx4","Hoxa11",#Allantois
  "Pcdh19","Dll1","Tceal6","Tbx6", #Somitic
  "Nkx1-2","Cdx4","Cdx2","Epha5",#NMP
  "Bhlha9",
  "Tfap2a","Grhl2","Grhl3",#Surface
  "Krt8","Krt18","Apela","Krtap17−1",#Gut
  "Rfx4" ,"Nrcam","Sfrp1",#Spinal cord
  "Sox2","Foxg1",
  
  "Isl1", "Meis1","Lefty2",##Pharyngeal mesoderm

  "Otx1","En1","Wnt8b","Six3",#Forebrain
  "Ptn",  "Sox2",  "Gbx2","Dmbx1", #midhind "En1",
  "Tubb2b",#
  "Nkx6-1","Mafb"#spinal cord
  
)

#Convert genes to gene IDs

geneIDs = getGeneID(genes)[[2]]

#Get coordinates of genes

TSS_peaks = read.table(paste0(wd,"PhD_BPS53/data/revision1/peaks_celltypes_inTSS.bed"),header=F,sep="\t")
TSS_peaks_df = data.frame(
  coord = paste0(as.character(TSS_peaks$V1),"_",as.character(TSS_peaks$V2),"_",as.character(TSS_peaks$V3)),
  gene = TSS_peaks$V8
)

#Subset coordinates of selected genes

TSS_peaks_sel = TSS_peaks_df[TSS_peaks_df$gene%in%geneIDs,]

#Create gene table with geneName and geneID

geneTable = getGeneID(genes)[[1]]
rownames(geneTable) = as.character(geneTable$Gene.ID)

#Subset counts of selected TSS
counts2 = counts[as.character(TSS_peaks_sel$coord),]
rownames(counts2)=as.character(TSS_peaks_sel$gene)

#Compute number of cells per cell type

n_celltype = table(m$final_clusters)

#Set column names with annotated cell types

counts2 = counts2[,names(n_celltype)]

colnames(counts2) = clustAnn_ATAC[names(n_celltype)]

#Compute frequency of accessibility

counts_freq = t(apply(counts2,1,function(x){x/n_celltype}))
colnames(counts_freq) = colnames(counts2)

#Standardise values

heatCountsStd <- t(apply(counts_freq, 1, function(x) x/max(x)))  # standarise


geneNames = getGeneName(rownames(heatCountsStd))[[2]]


#Set all parameters for heatmap
heatmapRedYelBlue <- c("#D73027", "#F46D43", "#FDAE61", "#FEE090",
                       "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4")
palette <- colorRampPalette(rev(heatmapRedYelBlue))


heatmapVir = as.character(read.csv(paste0(wd,"bin_data/viridis_100.tab"))[,1])
palette <- colorRampPalette((heatmapVir))

#Compute optimal reordering for heatmap
source(paste0(wd,"bin/Leaf_reordering_mod.R"))
den <- leaf.reordering(heatCountsStd,reorder = "both",corMethod="pearson")


#plot heatmap
pdf(paste0(wd,"PhD_BPS53/plots/revision1/snATACseq_revision1_genes_annotation_markerConsensus.pdf"),
    width=15,height=30)
gplots::heatmap.2(heatCountsStd[den$gene.clust$labels[den$gene.clust$order],den$cell.clust$labels[den$cell.clust$order]],
                  trace="none", 
                  col=palette, 
                  Colv = F, Rowv = F,  
                  labRow = anSeq::getGeneName(den$gene.clust$labels[den$gene.clust$order])[[2]],
                  #RowSideColors = cols, 
                  dendrogram = "none", density.info = 'none',
                  key=T,
                  #labCol=colnames(heatCountsVal),
                  margins= c(35,18),cexRow = 1.8,cexCol=2.5)
dev.off()


