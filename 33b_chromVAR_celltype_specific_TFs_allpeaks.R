###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to find top cell type-specific TF motifs.
###############################################################################


wd <- "/path/to/directory/"

#Define colour palettes
heatmapRedYelBlue <- c(
                       "#D73027", "#F46D43", "#FDAE61", "#FEE090",
                       "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4"
                       
                       )

palette <- colorRampPalette(rev(heatmapRedYelBlue))

all_colours = c(
  "Erythroid1" =  "#C72228",                              
  "Erythroid2" =  "#EF4E22", 
  "Erythroid3" =  "#f77b59",
  
  
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
all_colours_dot = all_colours
names(all_colours_dot) = gsub(" ",".",names(all_colours_dot))
names(all_colours_dot) = gsub("/",".",names(all_colours_dot))

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
####====================================
## ChromVAR
####====================================
#Read chromVAR results
meta = read.table(paste0(wd,"PhD_BPS53/data/revision1/snATACseq_embryo_revision1_metadata_afterAllQC_passQC_UMAP_cluster.txt"),
                  sep="\t",header=T)
rownames(meta)=meta$barcode
meta$final_clusters = paste0("cluster_",as.character(meta$final_clusters))


#read results from wilcoxon test
markers <- read.csv(paste0(wd,"PhD_BPS53/data/revision1/16_chromVAR/snATACseq_embryo_revision01_markerRegions_TFchromVar_allPeaks.csv"),
                    sep=",")
adjPval <- read.csv(paste0(wd,"PhD_BPS53/data/revision1/16_chromVAR/snATACseq_embryoTogether_markerRegions_pvalAdj_allPeaks.csv"),
                    sep=",")
logFC <- read.csv(paste0(wd,"PhD_BPS53/data/revision1/16_chromVAR/snATACseq_embryo_revision01_markerRegions_logFC_TFchromVar_allPeaks.csv"),
                  sep=",")
scores <- read.csv(paste0(wd,"PhD_BPS53/data/revision1/16_chromVAR/snATACseq_embryo_revision01_markerRegions_scores_TFchromVar_allPeaks.csv"),
                  sep=",")

colnames(scores) = clustAnn_ATAC[as.character(gsub("cluster_","",colnames(scores)))]
colnames(logFC) = clustAnn_ATAC[as.character(gsub("cluster_","",colnames(logFC)))]
colnames(adjPval) = clustAnn_ATAC[as.character(gsub("cluster_","",colnames(adjPval)))]
colnames(markers) = clustAnn_ATAC[as.character(gsub("cluster_","",colnames(markers)))]

chromVAR_Zscore = Matrix::readMM(paste0(wd,"PhD_BPS53/data/revision1/16_chromVAR/CHROMVAR_embryo_allPeaks_jointTogether_chromVAR_devScores_pwm_chromVar_sparse_raw_round4_peaksAll.mtx"))
colnames(chromVAR_Zscore) <- read.table(paste0(wd,"PhD_BPS53/data/revision1/16_chromVAR/CHROMVAR_embryo_allPeaks_jointTogether_chromVAR_devZscores_colnames_cellNames_pwm_chromVar_raw_round4_peaksAll.tab"))[,1]
rownames(chromVAR_Zscore) <- read.table(paste0(wd,"PhD_BPS53/data/revision1/16_chromVAR/CHROMVAR_embryo_allPeaks_jointTogether_chromVAR_devScores_rownames_motifs_pwm_chromVar_raw_round4_peaksAll.tab"))[,1]



####=======================================================
# Plot enrichment scores for GATA motifs in the UMAP
####=======================================================
GATA_motifs = rownames(chromVAR_Zscore)[grep("GATA",rownames(chromVAR_Zscore),ignore.case = T)]
chromVAR_GATA = chromVAR_Zscore[GATA_motifs,]

for (i in 1:nrow(chromVAR_GATA)){
  factor = rownames(chromVAR_GATA)[i]
  df = data.frame(x = meta[colnames(chromVAR_GATA),"umap_X"],
                  y=meta[colnames(chromVAR_GATA),"umap_Y"],
                  value=chromVAR_GATA[i,])
  df = df[order(df$value,decreasing=F),]
  
  
  
  tiff(paste0(wd,"PhD_BPS53/plots/revision1/UMAP_revision01_factor_",factor,".tiff"),
       width=7,height=7,res=300,units="in")
  
  p = ggplot(df, aes(x=x, y=y))
  plot = p + geom_point(aes(col = value))+
    theme_void() +
    
    scale_fill_gradient2(limits = c(-11,25),
                         high = "#de1404", mid = "#BFBFBF",
                           low = "#0664cf", midpoint = 0, space = "Lab",
                           na.value = "white", 
                           guide = "colourbar", aesthetics = "colour") +
    theme(legend.title = element_blank(),legend.position = "none") + 
    ggtitle("")
  print(plot)
  dev.off()
}



###Obtain top TFs

numbers = c()
TFs<-c()
TFmodules = list()
for (x in 1:ncol(markers)){
  m0 = as.character(markers[,x])
  p = adjPval[,x]
  l = logFC[,x]
  m = m0[p<0.0001]
  s = scores[,x][p<0.0001]

  
  TF= m[order(s,decreasing=T)]
  
  TF = TF[grep("ENSMUSG",TF)]
  TF= TF[1:15]
  
  numbers=c(numbers,length(TF))
  TFs <- c(TF,TFs)
}
names(numbers) = colnames(markers)


TFs <- unique(TFs)
length(TFs)


chromVAR_sel = chromVAR_Zscore[rownames(chromVAR_Zscore)%in% as.character(TFs),]
c <- as.character(meta[colnames(chromVAR_sel),"ann"])
heatCountsMean0 <- t(apply(chromVAR_sel,1,function(x){
  
  y <- aggregate(x,by=list(c),FUN=mean)
  y2 <- y$x[order(as.character(y$Group.1))]
  names(y2) <- as.character(y$Group.1)[order(as.character(y$Group.1))]
  
  
  return(y2)
}))





library(gplots)
###Tidy up TF names

d = strsplit(rownames(heatCountsMean0),"_")
TF_names=c()
for (i in d){
  if (substr(i[1],1,3) =="ENS" ){
    TF_names = c(TF_names,i[3])
  } else {
    TF_names = c(TF_names,i[4])
  }
}
type= c()
for (i in d){
  type = c(type,i[1])
}
discard = which(type %in% c("NP","Q8K439","Q8VHG7","XP"))


###Reformat values for heatmap plot
heatCountsMean = heatCountsMean0
rownames(heatCountsMean) = TF_names

heatCountsVal <- heatCountsMean[apply(heatCountsMean,1,function(x){var(x)>0}),]


#Perform optimal leaf reordering
source(paste0(wd,"bin_scripts/Leaf_reordering_mod.R"))
den <- leaf.reordering(heatCountsVal,reorder = "both",corMethod="pearson")

heatCountsStd <- t(apply(heatCountsVal, 1, function(x) x/max(x)))  # standarise

heatmapRedYelBlue <- c("#D73027", "#F46D43", "#FDAE61", "#FEE090",
                       "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4")
den$cell.clust$labels[den$cell.clust$order]
palette <- colorRampPalette(rev(heatmapRedYelBlue))

pdf(paste0(wd,"PhD_BPS53/plots/revision1/16_chromVAR/snATACseq_allTogether_chromVAR_std_peaksAll.pdf"),
    width=20,height=60)
gplots::heatmap.2((heatCountsStd[den$gene.clust$labels[den$gene.clust$order],den$cell.clust$labels[den$cell.clust$order]]), trace="none", 
                  col=palette, 
                  Colv = F, Rowv = F,  
                  #labRow = TF_names,
                  #RowSideColors = cols, 
                  dendrogram = "none", density.info = 'none',
                  key=T,breaks=seq(-1,1,0.1),
                  #labCol=colnames(heatCountsVal),
                  margins= c(35,18),cexRow = 1.8,cexCol=2.5)
dev.off()

pdf(paste0(wd,"PhD_BPS53/plots/revision1/16_chromVAR/snATACseq_allTogether_chromVAR_std_peaksAll_squeezed.pdf"),
    width=15,height=25)
gplots::heatmap.2((heatCountsStd[den$gene.clust$labels[den$gene.clust$order],den$cell.clust$labels[den$cell.clust$order]]), trace="none", 
                  col=palette, 
                  Colv = F, Rowv = F,  
                  #labRow = TF_names,
                  #RowSideColors = cols, 
                  dendrogram = "none", density.info = 'none',
                  key=T,breaks=seq(-1,1,0.1),
                  #labCol=colnames(heatCountsVal),
                  margins= c(35,18),cexRow = 1.5,cexCol=2.5)
dev.off()

heatCountsStdATAC = heatCountsStd
heatCountsMeanATAC = heatCountsMean

TFenriched = rownames(heatCountsStd)

#Save transcription factors enriched. This will be used in the following script for scRNA-seq
save(TFenriched,den,file=paste0(wd,"PhD_BPS53/data/revision1/16_chromVAR/snATACseq_enriched_TFs_top15_wilcoxon_allPeaks.rda"))

##==================================
## Compare to RNA
##  Note: Before this step, please execute 33c_scRNAseq_E825_analysis_TFexpression.R
##==================================


load(file=paste0(wd,"PhD_BPS53/data/revision1/16_chromVAR/ExpressionTFs_RNA_celltypes_wilcoxon_20190310_peaksAll.rda"))

#Homogeneise nomenclature
heatCountsMeanATAC =heatCountsMeanATAC[rownames(heatCountsMeanRNA),]
heatCountsMeanATAC = as.data.frame(heatCountsMeanATAC)
heatCountsMeanATAC$Erythroid1 =heatCountsMeanATAC[colnames(heatCountsMeanATAC)=="Erythroid"]
heatCountsMeanATAC$Erythroid2 =heatCountsMeanATAC[colnames(heatCountsMeanATAC)=="Erythroid"]
heatCountsMeanATAC$Erythroid3 =heatCountsMeanATAC[colnames(heatCountsMeanATAC)=="Erythroid"]

heatCountsMeanATAC[,"Intermediate mesoderm"] =heatCountsMeanATAC[colnames(heatCountsMeanATAC)=="Mixed mesoderm"]
heatCountsMeanATAC[,"ExE mesoderm"] =heatCountsMeanATAC[colnames(heatCountsMeanATAC)=="Mixed mesoderm"]

heatCountsMeanRNA = as.data.frame(heatCountsMeanRNA)
heatCountsMeanRNA[,"Forebrain"] =heatCountsMeanRNA[colnames(heatCountsMeanRNA)=="Forebrain/Midbrain/Hindbrain"]
heatCountsMeanRNA[,"Mid/Hindbrain"] =heatCountsMeanRNA[colnames(heatCountsMeanRNA)=="Forebrain/Midbrain/Hindbrain"]

heatCountsMeanATAC0 = heatCountsMeanATAC

heatCountsMeanATAC = heatCountsMeanATAC0[,(colnames(heatCountsMeanATAC0) %in% c("Intermediate mesoderm","ExE mesoderm","Forebrain/Midbrain/Hindbrain"))==F]
heatCountsMeanRNA0 = heatCountsMeanRNA

heatCountsMeanRNA = heatCountsMeanRNA0[,(colnames(heatCountsMeanRNA0) %in% c("Intermediate mesoderm","ExE mesoderm","Forebrain/Midbrain/Hindbrain"))==F]

celltypesOv = intersect(colnames(heatCountsMeanATAC),colnames(heatCountsMeanRNA))
heatCountsMeanRNA = heatCountsMeanRNA[,celltypesOv]
heatCountsMeanATAC = heatCountsMeanATAC[,celltypesOv]
heatCountsMeanATAC = as.matrix(heatCountsMeanATAC)
#pdf(paste0(wd,"PhD_BPS53/plots/Correlation_RNA_chromatin.pdf"),width=15,height=15)
par(mar=c(4,4,4,4),mfrow=c(8,10))


maxATAC=max(heatCountsMeanATAC)
minATAC=min(heatCountsMeanATAC)
maxRNA=max(heatCountsMeanRNA)
minRNA=min(heatCountsMeanRNA)

for (gene in rownames(heatCountsMeanRNA)){
  pdf(paste0(wd,"PhD_BPS53/plots/revision1/16_chromVAR/correlation_TF_RNA_chromatin_markerPeaks/Correlation_RNA_chromatin_",gene,".pdf"),
      width=5,height=5)
  par(mar=c(7,7,4,4))

  plot(unlist(unname(heatCountsMeanRNA[gene,])),unlist(unname(heatCountsMeanATAC[gene,])),
       pch=20,cex=2.5,cex.lab=2,cex.axis=1.5,cex.main=2.5,#main=gene,
        col=all_colours[colnames(heatCountsMeanRNA)],
       #ylab="Chromatin Accessibility",
       #xlab="RNA expression",
       ylab="",
       xlab="",
       ylim=c(minATAC,maxATAC)#,xlim=c(0,maxRNA)
       )
  #abline(h=0,lty=4)
  if (sum(heatCountsMeanRNA[gene,])>0){
    #abline(lm(heatCountsMeanATAC[gene,] ~ heatCountsMeanRNA[gene,]),col="red")
    
  }
  dev.off()
}

for (gene in rownames(heatCountsMeanRNA)){
  tiff(paste0(wd,"PhD_BPS53/plots/revision1/16_chromVAR/correlation_TF_RNA_chromatin_markerPeaks/Correlation_RNA_chromatin_",gene,".tiff"),
      width=5,height=5,res=300,units = "in")
  par(mar=c(7,7,4,4))
  
  plot(unlist(unname(heatCountsMeanRNA[gene,])),unlist(unname(heatCountsMeanATAC[gene,])),
       pch=20,cex=2.5,cex.lab=2,cex.axis=1.5,cex.main=2.5,#main=gene,
       col=all_colours[colnames(heatCountsMeanRNA)],
       #ylab="Chromatin Accessibility",
       #xlab="RNA expression",
       ylab="",
       xlab="",
       ylim=c(minATAC,maxATAC)#,xlim=c(0,maxRNA)
  )
  #abline(h=0,lty=4)
  if (sum(heatCountsMeanRNA[gene,])>0){
    #abline(lm(heatCountsMeanATAC[gene,] ~ heatCountsMeanRNA[gene,]),col="red")
    
  }
  dev.off()
}

pdf(paste0(wd,"PhD_BPS53/plots/revision1/16_chromVAR/TF_colorBar.pdf"),
    width=13,height=3)
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='',ylim=c(0,0.5),xlim=c(0,2))

xl <- 0.2; yb <- 0; xr <- 1.5; yt <- 0.2
rect(xleft=head(seq(xl,xr,(xr-xl)/50),-1), ybottom = yb, 
     xright = tail(seq(xl,xr,(xr-xl)/50),-1), ytop = yt,
     col=palette(50), border=palette(50), lwd=0.5, xpd=NA)
dev.off()#dev.off()

