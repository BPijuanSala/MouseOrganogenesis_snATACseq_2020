###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to define cell type-specific genomic regions and convey peak annotation
###############################################################################


library(foreach)
library(doMC)
registerDoMC(20)  #change the 2 to your number of CPU cores  
##--------------------------------------
#Computed in SERVER with SLURM
##--------------------------------------

wd<- "/path/to/directory/sample_pooled_preprocess_revision1/"
cat("\n reading meta\n")
meta = read.csv(paste0(wd,"14_visualisation/data/snATACseq_embryo_revision1_metadata_afterAllQC_passQC_UMAP_cluster.txt"),
                sep="\t")
rownames(meta)=meta$index

meta$final_clusters = paste0("cluster_",meta$final_clusters)



##===================
# Take counts per OCR
##===================
cat("\n reading OCR\n")

counts0 = read.table(paste0(wd,"14_visualisation/data/snATACseq_embryo_revision01_pooledNumber_OCR_perCelltype.csv"),
                     sep=",",header=T)
counts = counts0[,-1]
rownames(counts) = as.character(read.table(paste0(wd,"14_visualisation/data/snATACseq_embryo_revision01_pooledNumber_OCR_perCelltype_peakLabel.csv"),
                                           sep=",",header=T)[,2])
colnames(counts) = paste0("cluster_",as.character(read.table(paste0(wd,"14_visualisation/data/snATACseq_embryo_revision01_pooledNumber_OCR_perCelltype_celltypeLabel.csv"),
                                           sep=",",header=T)[,2]))




n_celltype = table(meta$final_clusters)
#names(n_celltype)[names(n_celltype)=="Mixed mesoderm"] = "ExE mesoderm"
##===================
# Compute fishers test one vs rest
##===================
cat("\n making functions\n")



evaluate = function(i,x,z){
  fishers_pval = matrix(0L,nrow=length(x),ncol=length(z))
  rownames(fishers_pval) = x
  colnames(fishers_pval) = z
  for (j in x){

    present_type = counts[j,i]
    absent_type = n_celltype[i]-present_type
    present_rest = sum(counts[j,colnames(counts)!=i])
    absent_rest =  sum(n_celltype[names(n_celltype)!=i])-present_rest
    matFish = matrix(0L,2,2)
    matFish[2,2]  =present_type
    matFish[2,1]  =absent_type
    matFish[1,2]  =present_rest
    matFish[1,1]  =absent_rest
    #With alternative="greater" you're asking, "Is there more signal in the upper right and lower left corners than expected?"
    
    fishers_pval[j,i] = fisher.test(matFish,alternative = "greater")$p.value
    
    
  }
  return(fishers_pval[,i])
}

cat("\n create matrix combine\n")

matrix<-foreach(i=colnames(counts),.combine=rbind) %dopar% {
  cat(i)
  evaluate(i,x=rownames(counts),z=colnames(counts))
}
fishers_pval = matrix
save(matrix,file=paste0(wd,"14_visualisation/data/snATACseq_revision1_peaks_enriched_fishersTest.rda"))



##===================
# Compute fishers test pairwise
##===================

load(file=paste0(wd,"14_visualisation/data/snATACseq_revision1_peaks_enriched_fishersTest.rda"))
fishers_pval = t(matrix)

colnames(fishers_pval) = colnames(counts)
rownames(fishers_pval) =rownames(counts)

n_celltype = table(meta$final_clusters)


evaluate = function(x,testSites,others,counts,n_celltype,celltype,fishers_test){
  
  for (i in 1:length(testSites)){
    cat(i,"...")
    
    for (j in 1:length(others)){
      
      otherCelltype = others[j]
      present_type = counts[testSites[i],celltype]
      absent_type = n_celltype[celltype]-present_type
      present_rest = sum(counts[testSites[i],otherCelltype])
      absent_rest =  sum(n_celltype[otherCelltype])-present_rest
      matFish = matrix(0L,2,2)
      matFish[2,2]  =present_type
      matFish[2,1]  =absent_type
      matFish[1,2]  =present_rest
      matFish[1,1]  =absent_rest
      #With alternative="greater" you're asking, "Is there more signal in the upper right and lower left corners than expected?"
      
      fishers_test[testSites[i],otherCelltype] = fisher.test(matFish,alternative = "greater")$p.value
      
    }
    
  }
  return(fishers_test)
  
}

fishers_Celltype<-foreach(x=1:ncol(fishers_pval)) %dopar% {
  celltype = colnames(fishers_pval)[x]
  cat(celltype,"\n")
  others = colnames(fishers_pval)[colnames(fishers_pval) != celltype]
  st = p.adjust(fishers_pval[,x],method="bonferroni")
  testSites = names(st)[st < 0.01]
  
  fishers_test = matrix(0L,nrow=length(testSites),ncol=length(others))
  rownames(fishers_test) = testSites
  colnames(fishers_test) = others
  
  evaluate(x,testSites,others,counts,n_celltype,celltype,fishers_test)
  
  
} 


save(fishers_Celltype,file=paste0(wd,"14_visualisation/data/peaks_enriched_fishersTest_pairwise_selected_round4.rda"))



##--------------------------------------
#Computed LOCALLY
##--------------------------------------
library(ggplot2)
wd <- "/path/to/directory/"
#================
# Clusters and colours
#================


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

all_colours = c(
  "Allantois" = "#532C8A",
  "Cardiomyocytes" =  "#B51D8D",  
  "Mixed mesoderm" =  "#ab80b7",# E8.25 mixed mesoderm   
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
names(all_colours_dot) =  gsub(" ",'.',names(all_colours))
names(all_colours_dot) =  gsub('/','.',names(all_colours_dot))
###################################
# Define marker regions
###################################

# Read counts per cluster
counts0 = read.table(paste0(wd,"PhD_BPS53/data/revision1/snATACseq_embryo_revision01_pooledNumber_OCR_perCelltype.csv"),
                     sep=",",header=T)
counts = counts0[,-1]
rownames(counts) = as.character(read.table(paste0(wd,"PhD_BPS53/data/revision1/snATACseq_embryo_revision01_pooledNumber_OCR_perCelltype_peakLabel.csv"),
                                           sep=",",header=T)[,2])
colnames(counts) = paste0("cluster_",as.character(read.table(paste0(wd,"PhD_BPS53/data/revision1/snATACseq_embryo_revision01_pooledNumber_OCR_perCelltype_celltypeLabel.csv"),
                                           sep=",",header=T)[,2]))


meta = read.csv(paste0(wd,"PhD_BPS53/data/revision1/snATACseq_embryo_revision1_metadata_afterAllQC_passQC_UMAP_cluster.txt"),
                sep="\t")
rownames(meta)=meta$barcode
meta$final_clusters = paste0("cluster_",meta$final_clusters)



n_celltype = table(meta$final_clusters)

#load fisher's exact test results
load(file=paste0(wd,"PhD_BPS53/data/revision1/14_visualisation/snATACseq_revision1_peaks_enriched_fishersTest.rda"))

fishers_pval = t(matrix)
colnames(fishers_pval)=colnames(counts)


load(file=paste0(wd,"PhD_BPS53/data/revision1/14_visualisation/peaks_enriched_fishersTest_pairwise_selected_round4.rda"))
names(fishers_Celltype)=colnames(fishers_pval)


#Get marker regions setting by up a threshold
fishers_Celltype_qval = fishers_Celltype
celltype_markers = list()


for (i in 1:length(fishers_Celltype_qval)){
  cat(names(fishers_Celltype)[i],"...")
  fishers_qval = fishers_Celltype[[i]]
  for (x in 1:nrow(fishers_qval)){
    fishers_qval[x,] = p.adjust(fishers_qval[x,],method="bonferroni")
    
  }
  fishers_Celltype_qval[[i]] = fishers_qval
  celltype_markers[[names(fishers_Celltype)[i]]] = rownames(fishers_qval)[apply(fishers_qval,1,function(x){length(which(x>=1e-10))<=9})]
  dat = apply(counts[sample(celltype_markers[[names(fishers_Celltype)[i]]],10,replace=T),],1,function(x){x/n_celltype[colnames(counts)]})
  barplot(dat,beside = T,col=all_colours[clustAnn_ATAC[gsub("cluster_","",rownames(dat))]],las=1,main=paste(names(fishers_Celltype)[i],length(celltype_markers[[names(fishers_Celltype)[i]]])))
  
}
lengths(celltype_markers)


#########====================

#Accessibility in genomic regions
access = read.csv(paste0(wd,"PhD_BPS53/data/revision1/14_visualisation/snATACseq_embryo_revision1_cellsAsVars_TFIDF_with_cellsObs_accessibility.csv"))

umapPeakObs = read.csv(paste0(wd,"PhD_BPS53/data/revision1/14_visualisation/snATACseq_embryo_revision1_cellsAsVars_TFIDF_with_cellsObs_UMAP.csv"))

metaPeak = data.frame(
  "peak" = access[,1],
  "accessibility"=access$accessibility,
  "accessibility_log"=log(access$accessibility+1),
  "accessibility_ratio"=access$accessibility/19453,
  "umapX"=umapPeakObs$X0,
  "umapY"=umapPeakObs$X1
)
rownames(metaPeak)=metaPeak$peak




# Plot marker regions in UMAP
for (i in 1:length(celltype_markers)){
  c = names(celltype_markers)[i]
  c = gsub(' ',"_",c)
  c = gsub('/',"_",c)
  
  tiff(paste0(wd,"PhD_BPS53/plots/revision1/20190309_all_together_annotatedPeaks_classification_markersClusters_stringentQC_refined_fishersTest_pairwise",c,".tiff"),
      width=5,height=5.2,res=300,units="in")

  plot(metaPeak$umapX,metaPeak$umapY,pch=20,col="#dbd9d9",
       cex=0.1,#col="black",
       xaxt="n",yaxt="n",ylab="",xlab="",bty="n"#,main=names(celltype_markers)[i]
       )
  idx = rownames(metaPeak) %in% as.character(celltype_markers[[i]])
  points(metaPeak[idx,"umapX"],metaPeak[idx,"umapY"],
         pch=20,col=unname(all_colours[clustAnn_ATAC[as.character(gsub("cluster_","",c))]]),cex=0.25)
  dev.off()
}
#names(celltype_markers) = clustAnn_ATAC[gsub("cluster_","",names(celltype_markers))]


# Define ubiquitous regions
hist(log(metaPeak$accessibility_ratio),breaks=100)
abline(v=log(0.25))

ubiq = metaPeak$accessibility_ratio>0.25


tiff(paste0(wd,"PhD_BPS53/plots/revision1/20190309_all_together_annotatedPeaks_classification_markersClusters_stringentQC_refined_fishersTest_pairwise_ubiquitous_10lineages.tiff"),
     width=5,height=5.2,res=300,units="in")
plot(metaPeak$umapX,metaPeak$umapY,pch=20,col="#dbd9d9",cex=0.1,
     xaxt="n",yaxt="n",ylab="",xlab="",bty="n",main="")
points(metaPeak[ubiq,"umapX"],metaPeak[ubiq,"umapY"],
       pch=20,col="black",cex=0.25)
dev.off()


# Plot accessibility
cols = c("#0F0F81","#21AB77","#FBFB95","#FFA500","#B01414")
mixRamp <- colorRampPalette(cols)
metaPeak$accessibility_ratio_log = log(metaPeak$accessibility_ratio)
pdf(paste0(wd,"PhD_BPS53/plots/revision1/20190502_all_together_peaks_UMAP_accessibility.pdf"),
    width=10,height=7)
p = ggplot(metaPeak, aes(x=umapX, y=umapY))
p + geom_point(aes(col = accessibility_ratio_log),size=0.01)+
  theme_void() +
  scale_colour_gradientn(colours = cols) #+ ggtitle("Read depth")

dev.off()

pdf(paste0(wd,"PhD_BPS53/plots/revision1/accessibility_colorBar.pdf"),
    width=13,height=3)
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='',ylim=c(0,0.5),xlim=c(0,2))

xl <- 0.2; yb <- 0; xr <- 1.5; yt <- 0.2
rect(xleft=head(seq(xl,xr,(xr-xl)/50),-1), ybottom = yb, 
     xright = tail(seq(xl,xr,(xr-xl)/50),-1), ytop = yt,
     col=mixRamp(50), border=mixRamp(50), lwd=0.5, xpd=NA)
dev.off()

tiff(paste0(wd,"PhD_BPS53/plots/revision1/20190502_all_together_peaks_UMAP_accessibility.tiff"),
    width=7,height=7,res=300,units="in")
p = ggplot(metaPeak, aes(x=umapX, y=umapY))
p + geom_point(aes(col = accessibility_ratio_log),size=0.01)+
  theme_void() + theme(legend.title = element_blank()) +
  scale_colour_gradientn(colours = cols) #+ ggtitle("Read depth")

dev.off()

range(metaPeak$accessibility_ratio)

############GET CELL TYPE MARKERS
for (i in 1:length(celltype_markers)){
  clust = names(celltype_markers)[i]
  coords = celltype_markers[[i]]
  write.table(coords,file=paste0(wd,"PhD_BPS53/data/revision1/10_peaks_all/celltypes_coords/",clust,".txt"),
              sep="\n",quote=F,row.names=F,col.names=F)

  }

#Write files

markers = sort(unique(unname(unlist(celltype_markers))))
#write.table(markers,sep="\n",quote=F,col.names=F,row.names=F,
#            file = paste0(wd,"PhD_BPS53/data/revision1/14_visualisation/snATACseq_peaks_markerRegions.txt"))
library(reshape2)

markers = melt(celltype_markers)
markers_df = data.frame(
  "chr"=gsub("_.*","",markers$value),
  "chr"=gsub("_.*","",gsub("chr[0-9]*_","",markers$value)),
  "end"=gsub(".*_","",markers$value),
  "celltype" = markers$L1
)
#write.table(markers_df,sep="\t",quote=F,col.names=F,row.names=F,
 #          file = paste0(wd,"PhD_BPS53/data/revision1/14_visualisation/snATACseq_peaks_markerRegions.txt"))






##############Explore connections as network/heatmap based on shared marker peaks
markers = sort(unique(unname(unlist(celltype_markers))))


# Make matrix with cell type-specific peaks and cell types
markersMat = matrix(0L,nrow=length(markers),ncol=length(celltype_markers))
colnames(markersMat) =names(celltype_markers)
rownames(markersMat) = markers

for (i in rownames(markersMat)){
  markersMat[i,] = as.numeric(unlist(lapply(celltype_markers,function(x){i %in% x})))
}



colnames(markersMat) = unname(clustAnn_ATAC[as.character(gsub("cluster_","",colnames(markersMat)))])
tab = table(rowSums(markersMat))

pdf(paste0(wd,"PhD_BPS53/plots/revision1/Celltype_specific_shared.pdf"),
    width=4,height=4)
barplot(log(tab),ylim=c(0,10),col="#0869a1")

dev.off()

#How much each lineage shares
sharing = matrix(0L,nrow=length(celltype_markers),ncol=6)
rownames(sharing) =colnames(markersMat)
colnames(sharing) = 1:6
for (i in 1:nrow(sharing)){
  tmp = markersMat[markersMat[,i]>0,]
  tab = table(rowSums(tmp))
  sharing[i,names(tab)] = tab
}


# Compute adjacency

markersMat3 = markersMat
range(colSums(markersMat3))
df1.mt <- as.matrix(markersMat3)
df1.adj <- t(df1.mt) %*% df1.mt
df2.adj = df1.adj
pairs = t(combn(colnames(markersMat3),2))

# Compute union

df2.union = df2.adj

for (i in 1:nrow(pairs)){
  pair1 = pairs[i,1]
  pair2 = pairs[i,2]
  df2.union[pair1,pair2] = df2.union[pair2,pair1] = length(which(rowSums(markersMat3[,c(pair1,pair2)])>0))
  
}

# Compute Jaccard similarity index
df2.adj.norm = (df2.adj  / df2.union)


# Set diagonal to 0

diag(df2.adj.norm) <- 0


# Create heatmap

heatmapRedYelBlue <- c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090",
                       "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")

palette <- colorRampPalette(rev(heatmapRedYelBlue))


pdf(paste0(wd,"PhD_BPS53/plots/revision1/heatmap_shared_celltype_specific_peaks.pdf"),
    width=7,height=7)
gplots::heatmap.2(df2.adj.norm, trace="none", 
                  col=palette, 
                  Colv = T,
                  Rowv = T, 
                  #ColSideColors = c, 
                  #RowSideColors = cols, 
                  dendrogram = "none", density.info = 'none',
                  #labRow = g,
                  key=T#labCol=""
)
dev.off()



### If doing a network...

# Set a threshold
hist(df2.adj.norm,col="black",ylim=c(0,50),breaks=70)
abline(v=0.03,col="red")


df2.adj.norm[df2.adj.norm < 0.03] = 0

# Create graph
library(igraph)
g <- graph_from_adjacency_matrix(df2.adj.norm,weighted=TRUE,diag=F,mode="undirected")
plot(g)

# Export graph
gexf = rgexf::igraph.to.gexf(g)
rgexf::write.gexf(nodes = gexf$nodes[,c(1,2)], 
                  edges = gexf$edges[,c(2,3)], 
                  edgesWeight = gexf$edges[,c("weight")],
                  nodesAtt = (as.data.frame(unique(gexf$nodes$label))), 
                  output = paste0(wd,"PhD_BPS53/data/revision1/14_visualisation/shared_edges_norm_1min10_9lineages_pairwise.gexf"))





#############peak Annotation

peakAnn = read.csv(paste0(wd,"PhD_BPS53/data/revision1/10_peaks_all/peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_annotated_assignGenomeAnnotation_withPeaks.bed"),
                   sep="\t",header=F)
colnames(peakAnn)=c("Chr","Start","End","ID","ID","noIdea","Annotation")
rownames(peakAnn)<-paste0("chr",peakAnn$Chr,"_",(peakAnn$Start),"_",peakAnn$End)
peakAnn =peakAnn[order(peakAnn$Chr,peakAnn$Start),]
as.character(unique(peakAnn$Annotation))
peakAnn$Annotation.General = gsub("\\(.*$","",gsub("--.*$","",gsub(" .*$","",peakAnn$Annotation)))
peakAnn$Annotation.General[peakAnn$Annotation.General =="promo"] = "promoter-TSS"
unique(peakAnn$Annotation.General)#TTs: transcription termination sites

peakAnn$ID = rownames(peakAnn)

annCol = c(
  "Intergenic"   = "red",
  "intron" = "orange",
  "promoter-TSS" ="blue",
  "TSS" ="blue",
  
  "TTS"="gold",          
  "exon"  = "lightblue"#,       
)


names(celltype_markers) = clustAnn_ATAC[as.character(gsub("cluster_","",names(celltype_markers)))]
celltype_markers$Ubiquitous = rownames(metaPeak)[metaPeak$accessibility_ratio>0.25]

lengths(celltype_markers)

##======================================================================
# Annotate marker regions
##======================================================================

listReg=list()



for (i in 1:length(celltype_markers)){
  celltype = names(celltype_markers)[i]
  reg = as.character(celltype_markers[[i]])
  for (j in 1:length(reg)){
    r= reg[j]
    
    
    ann = as.character(peakAnn[r,"Annotation.General"])
    
    listReg[[paste0(celltype,"_",r)]] = c(r,ann,celltype)
    
  }
  
}

df = t(do.call(data.frame,listReg))
colnames(df)=c("ID","ann","celltype")
df <- as.data.frame(df)
table(df$ID)
df[df$celltype=="Endothelium",]
df2 = df
total = cbind(peakAnn[,c("ID","Annotation.General")],rep("Total",nrow(peakAnn)))
colnames(total)=c("ID","ann","celltype")

df2 = rbind(df2,total)
tab1 = table(df2$ann,df2$celltype)

tab2=apply(tab1,2,function(x){x/sum(x)})
tab2 = tab2[,c(  "Total",
                 "Erythroid", "Endothelium" ,"Allantois", "Mixed mesoderm",
                 "Mesenchyme","Pharyngeal mesoderm","Cardiomyocytes",
                 "Paraxial mesoderm" ,"Somitic mesoderm" ,"NMP" ,
                 "Spinal cord","Forebrain", "Mid/Hindbrain",
                 "Neural crest"   , "Surface ectoderm" , "Notochord", "Gut",
                 "ExE endoderm" , "Ubiquitous")]
pdf(paste0(wd,"PhD_BPS53/plots/revision1/Celltype_specific_peaks/20190520_all_together_annotatedPeaks_classification.pdf"),
    width=7,height=3.5)
par(mar=c(10,5,2,2))
barplot(tab2[,order(tab2[c("promoter-TSS"),],-tab2["Intergenic",])],col=annCol[rownames(tab2)],
        las=2,horiz=F,xlab="")
dev.off()




#***************************************************************************
##===========================
# ALTERNATIVE to code above to annotate marker regions more refined
# This uses a variable "peakGene" generated in the code below.
##===========================
load(file=paste0(wd,"PhD_BPS53/data/revision1/10_peaks_all/peakGene.rda"))
peakGene[as.character(peakGene$distance_from_TSS) =="TSS" ,"Annotation.General"] = "TSS"
peakGene[peakGene$accessibility_ratio>0.25  ,"celltype_specificity"] = "Ubiquitous"

peakGene_uniq = peakGene[!duplicated(peakGene[,c('peakID')]),]
rownames(peakGene_uniq)=peakGene_uniq$peakID

listReg=list()



for (i in 1:length(celltype_markers)){
  celltype = names(celltype_markers)[i]
  reg = as.character(celltype_markers[[i]])
  for (j in 1:length(reg)){
    r= reg[j]
    
    
    ann = as.character(peakGene_uniq[r,"Annotation.General"])
    
    listReg[[paste0(celltype,"_",r)]] = c(r,ann,celltype)
    
  }
  
}

df = t(do.call(data.frame,listReg))
colnames(df)=c("ID","ann","celltype")
df <- as.data.frame(df)
table(df$ID)
df2 = df
total = cbind(peakGene_uniq[,c("ID","Annotation.General")],rep("Total",nrow(peakGene_uniq)))
colnames(total)=c("ID","ann","celltype")

df2 = rbind(df2,total)
tab1 = table(df2$ann,df2$celltype)

tab2=apply(tab1,2,function(x){x/sum(x)})
tab2 = tab2[,c(  "Total",
                 "Erythroid", "Endothelium" ,"Allantois", "Mixed mesoderm",
                 "Mesenchyme","Pharyngeal mesoderm","Cardiomyocytes",
                 "Paraxial mesoderm" ,"Somitic mesoderm" ,"NMP" ,
                 "Spinal cord","Forebrain", "Mid/Hindbrain",
                 "Neural crest"   , "Surface ectoderm" , "Notochord", "Gut",
                 "ExE endoderm" , "Ubiquitous")]
pdf(paste0(wd,"PhD_BPS53/plots/revision1/Celltype_specific_peaks/20190520_all_together_annotatedPeaks_classification.pdf"),
    width=7,height=3.5)
par(mar=c(10,5,2,2))
barplot(tab2[,order(tab2[c("TSS"),],-tab2["Intergenic",])],col=annCol[rownames(tab2)],
        las=2,horiz=F,xlab="")
dev.off()

##***************************************************************************




metaPeak$annotation = peakAnn[rownames(metaPeak),"Annotation.General"]

#write.table(metaPeak,paste0(wd,"PhD_BPS53/data/revision1/10_peaks_all/metadata_genomicRegions_final.txt"),
#        quote=F,sep="\t",col.names=T,row.names=F)





##############################################################
# Genes assigned
##############################################################

peakGene = read.csv(paste0(wd,"PhD_BPS53/data/revision1/10_peaks_all/peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist_nochr_geneAssignment_distance_plusUnmapped.bed"),
                    sep="\t",header=F)
colnames(peakGene) = c("peak_chr","peak_start","peak_end","geneID","geneName","strand","distance_from_TSS")

peakGene$peakID = paste0("chr",peakGene$peak_chr,"_",peakGene$peak_start,"_",peakGene$peak_end)
#TSS (by default defined from -1kb to +100bp)
#TTS (by default defined from -100 bp to +1kb)
#Distance is from TSS. If distance is TSS: TSS - 1kb
#If it's gene: It's stricly within the gene.
peakGene[,"Annotation.General"] = peakAnn$Annotation.General[match(peakGene$peakID,peakAnn$ID)]


peakGene[(peakGene$distance_from_TSS =="gene" & peakGene$Annotation.General=="promoter-TSS"),"distance_from_TSS"] = "TSS"
peakGene[(peakGene$Annotation.General =="promoter-TSS"),"Annotation.General"] = "TSS"
peakGene[(peakGene$distance_from_TSS =="TSS"),"Annotation.General"] = "TSS"

peakGene[unique(peakGene$peakID),]








table(peakGene$Annotation.General)

markers_m = reshape::melt(celltype_markers)
colnames(markers_m) = c("peakID","celltype")

peakGene[,"celltype_specificity"] = rep(NA,nrow(peakGene))
for (i in 1:length(unique(markers_m$celltype))){
  top=unique(markers_m$celltype)[i]
  coord = as.character(markers_m$peakID)[as.character(markers_m$celltype)==top]
  #idx = unlist(lapply(coord, function(x) {which(as.character(peakGene$peakID)==x)}))
  idx = which(!is.na(match(as.character(peakGene$peakID), coord)))
  peakGene[idx,"celltype_specificity"] = paste0(peakGene[idx,"celltype_specificity"],";",top,"")
  
}
peakGene[peakGene$peakID =="chr4_115075466_115075967",]
peakGene$celltype_specificity = gsub("NA;","",peakGene$celltype_specificity)

table(peakGene$celltype_specificity)

peakGene[,"topic"] = rep(NA,nrow(peakGene))
for (i in 1:100){
  top=i
  tab = rownames(read.table(paste0(wd,"PhD_BPS53/data/revision1/13_cisTopic/topics_regions/cisTopicObject_50_100_afterReadsInPeaksQC_24_Topics_Topic",top,".txt"),
                            sep="\t",header=T))
  coord = gsub("-","_",gsub(":","_",tab))
  idx = which(!is.na(match(as.character(peakGene$peakID), coord)))
  
  peakGene[idx,"topic"] = paste0(peakGene[idx,"topic"],";","Topic",top,"")
  
  
  
}
peakGene$topic = gsub("NA;","",peakGene$topic)
peakGene$topic_stringent = peakGene$topic
peakGene$topic_stringent[grep(";Topic",peakGene$topic_stringent)] = "Nonspecific"
peakGene = peakGene[,c("peakID",'peak_chr','peak_start','peak_end','Annotation.General',
                       'distance_from_TSS','geneName','geneID','strand','celltype_specificity',
                       'topic','topic_stringent')]


peakGene[,"accessibility"] = metaPeak$accessibility[match(peakGene$peakID,metaPeak$peak)]
peakGene[,"accessibility_log"] = metaPeak$accessibility_log[match(peakGene$peakID,metaPeak$peak)]
peakGene[,"accessibility_ratio"] = metaPeak$accessibility_ratio[match(peakGene$peakID,metaPeak$peak)]
peakGene[,"umap_X"] = metaPeak$umapX[match(peakGene$peakID,metaPeak$peak)]
peakGene[,"umap_Y"] = metaPeak$umapY[match(peakGene$peakID,metaPeak$peak)]

##===================================
# patterns endothelium - allantois
##===================================
load(file=paste0(wd,"PhD_BPS53/data/revision1/18_endothelium_analysis/cluster_DPT/snATACseq_endothelium_eryAl_DPT_clusters.rda"))
clust = clusters$clusters
peakGene[,"Pattern_endothelium"] = rep(NA,nrow(peakGene))
for (i in sort(unique(clusters$clusters))){
  coord = names(clust)[clust == i]
  idx = which(!is.na(match(as.character(peakGene$peakID), coord)))
  
  peakGene[idx,"Pattern_endothelium"] = i
  
}

#save(peakGene,file=paste0(wd,"PhD_BPS53/data/revision1/10_peaks_all/peakGene.rda"))

load(file=paste0(wd,"PhD_BPS53/data/revision1/10_peaks_all/peakGene.rda"))
peakGene[as.character(peakGene$distance_from_TSS) =="TSS" ,"Annotation.General"] = "TSS"
peakGene[peakGene$accessibility_ratio>0.25  ,"celltype_specificity"] = "Ubiquitous"
#save(file=paste0(wd,"PhD_BPS53/data/revision1/10_peaks_all/peakGene_refined.rda"))

#write.table(peakGene,paste0(wd,"PhD_BPS53/data/revision1/10_peaks_all/metadata_genomicRegions_final_peakGene.txt"),
#       quote=F,sep="\t",col.names=T,row.names=F)
paper = "/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/My_papers_submitted/20190519_PijuanSala_etal_snATACseq/02_Submission_NatCellBiol/02_revision_01/05_Tables/"

#write.table(peakGene,paste0(paper,"Supplementary_Table6_OCRs.txt"),
#       quote=F,sep="\t",col.names=T,row.names=F)

###################

##=========================
### SURFACE VS GUT
##=========================

int = intersect(as.character(celltype_markers$cluster_5),as.character(celltype_markers$cluster_7))
clust7only = setdiff(celltype_markers$cluster_7,int)
clust5only = setdiff(celltype_markers$cluster_5,int)
length(int)
length(clust7only)
length(clust5only)


library("org.Mm.eg.db")
library("GO.db")
library("topGO")
library("GOstats")

genesTest =unique(as.character(peakGene[match(int,peakGene$peakID),"geneID"]))

a <- sapply(genesTest, function(x) exists(x, org.Mm.egENSEMBL2EG))
my.ensmusg.existed <- genesTest[a]

xx <- as.list(org.Mm.egENSEMBL2EG)
EntrezGenes <- unlist(xx[my.ensmusg.existed])
UniverseEntrezGenes <- unlist(xx)
go_object <- as.list(org.Mm.egGO2EG)

params <- new('GOHyperGParams',
              geneIds=unique(EntrezGenes),
              universeGeneIds=unique(UniverseEntrezGenes),
              ontology='BP',
              pvalueCutoff=0.001,
              conditional=F,
              testDirection='over',
              annotation="org.Mm.eg.db"
)

hgOver <- hyperGTest(params)

result <- summary(hgOver)
result$fdr = p.adjust(result$Pvalue,method = "BH")
result = result[order(result$fdr,decreasing=F),]
head(result,20)#intersection


pdf(paste0(wd,"PhD_BPS53/plots/revision1/Celltype_specific_peaks/GO_terms_shared_surface_gut.pdf"),
    width=6,height=3)
par(mar=c(2,20,2,2))
barplot(-log(result$fdr[1:15]),main="", xlab="-log(FDR)", horiz=TRUE,
        names.arg=result$Term[1:15],las=1,col="black",xlim=c(0,15),cex.lab=2)
dev.off()


write.table(result,paste0(wd,"PhD_BPS53/data/revision1/GO_terms_gut_surface/intersection.txt"),
            sep="\t",quote=F,col.names=T,row.names=F)

#### CLUSTER 5 - Surface ectoderm

genesTest =unique(as.character(peakGene[match(clust5only,peakGene$peakID),"geneID"]))

a <- sapply(genesTest, function(x) exists(x, org.Mm.egENSEMBL2EG))
my.ensmusg.existed <- genesTest[a]

xx <- as.list(org.Mm.egENSEMBL2EG)
EntrezGenes <- unlist(xx[my.ensmusg.existed])
UniverseEntrezGenes <- unlist(xx)
go_object <- as.list(org.Mm.egGO2EG)

params <- new('GOHyperGParams',
              geneIds=unique(EntrezGenes),
              universeGeneIds=unique(UniverseEntrezGenes),
              ontology='BP',
              pvalueCutoff=0.001,
              conditional=F,
              testDirection='over',
              annotation="org.Mm.eg.db"
)

hgOver <- hyperGTest(params)

result_clust5 <- summary(hgOver)

result_clust5$fdr = p.adjust(result_clust5$Pvalue,method = "BH")
result_clust5 = result_clust5[order(result_clust5$fdr,decreasing=F),]
head(result_clust5,20)

pdf(paste0(wd,"PhD_BPS53/plots/revision1/Celltype_specific_peaks/GO_terms_surfaceEcto.pdf"),
    width=6,height=3)
par(mar=c(2,20,2,2))
barplot(-log(result_clust5$fdr[1:15]),main="", xlab="-log(FDR)", horiz=TRUE,
        names.arg=result_clust5$Term[1:15],las=1,col=all_colours["Surface ectoderm"],
        xlim=c(0,26),cex.lab=2)
dev.off()
write.table(result_clust5,paste0(wd,"PhD_BPS53/data/revision1/GO_terms_gut_surface/cluster5_GO.txt"),
            sep="\t",quote=F,col.names=T,row.names=F)


#### CLUSTER 7 - Gut
genesTest =unique(as.character(peakGene[match(clust7only,peakGene$peakID),"geneID"]))

a <- sapply(genesTest, function(x) exists(x, org.Mm.egENSEMBL2EG))
my.ensmusg.existed <- genesTest[a]

xx <- as.list(org.Mm.egENSEMBL2EG)
EntrezGenes <- unlist(xx[my.ensmusg.existed])
UniverseEntrezGenes <- unlist(xx)
go_object <- as.list(org.Mm.egGO2EG)

params <- new('GOHyperGParams',
              geneIds=unique(EntrezGenes),
              universeGeneIds=unique(UniverseEntrezGenes),
              ontology='BP',
              pvalueCutoff=0.001,
              conditional=F,
              testDirection='over',
              annotation="org.Mm.eg.db"
)

hgOver <- hyperGTest(params)

result_clust7 <- summary(hgOver)

result_clust7$fdr = p.adjust(result_clust7$Pvalue,method = "BH")
result_clust7 = result_clust7[order(result_clust7$fdr,decreasing=F),]
head(result_clust7,20)


write.table(result_clust7,paste0(wd,"PhD_BPS53/data/revision1/GO_terms_gut_surface/cluster7_GO.txt"),
            sep="\t",quote=F,col.names=T,row.names=F)


pdf(paste0(wd,"PhD_BPS53/plots/revision1/Celltype_specific_peaks/GO_terms_Gut.pdf"),
    width=6,height=3)
par(mar=c(2,20,2,2))
barplot(-log(result_clust7$fdr[1:15]),main="", xlab="-log(FDR)", horiz=TRUE,
        names.arg=result_clust7$Term[1:15],las=1,col=all_colours["Gut"],
        xlim=c(0,45),cex.lab=2)
dev.off()

head(result_clust7,20)






###############################################################
# GATA binding sites in topics

d = load(paste0(wd,"PhD_BPS53/data/revision1/16_chromVAR/CHROMVAR_motif_ix_pwm_chromVar_peaksAll.rda"))

mat = motif_ix@assays$data$motifMatches
rownames(mat) = rownames(motif_ix)
colnames(mat) = colnames(motif_ix)

GATA_TF = colnames(mat)[grep("Gata",colnames(mat))]
mat_GATA = mat[,GATA_TF]
GATA_reg = rownames(mat_GATA)[Matrix::rowSums(mat_GATA) > 0]


length(GATA_reg)
topics_GATA = sort(table(peakGene[match(GATA_reg,peakGene$peakID),"topic_stringent"]),decreasing=T)
hist(topics_GATA)
hist(topics_GATA[2:length(topics_GATA)],breaks=50)
peakGene[match(GATA_reg,peakGene$peakID),"topic_stringent"]

y=c()
for (i in 1:100){
  Topic1_motif_enr <- RcistargetOutput[[paste0("Topic",i)]]
  x = grep("Gata",Topic1_motif_enr[,"TF_highConf"])
  if (length(x)>0){
    y = c(y,paste0("Topic",i))
  }
  
}
names(topics_GATA)=gsub(";","",names(topics_GATA))


pdf(paste0(wd,"PhD_BPS53/plots/revision1/NumberGATA_per_topic.pdf"),
    width=5,height=4)
hist(sort(topics_GATA[y],decreasing=T),breaks=15,col="black",
     main="",xlab="Number of regions with GATA binding sites", ylab="Number of topics")
dev.off()


###################################################################
# TRANSGENICS - check peaks
###################################################################
load(file=paste0(wd,"PhD_BPS53/data/revision1/10_peaks_all/peakGene.rda"))

peakGene[as.character(peakGene$geneName)=="Tal1",]

peakGene[as.character(peakGene$geneName)=="Erg",]

#Erg cell type-specific peak
peakGene[as.character(peakGene$peakID)=="chr16_95438997_95439532",]

#Fli1 -15kb
peakGene[as.character(peakGene$geneName)=="Fli1",]
peakGene[as.character(peakGene$peakID)=="chr9_32556231_32556843",]

#Flt1 +67kb
peakGene[as.character(peakGene$geneName)=="Flt1",]
peakGene[as.character(peakGene$peakID)=="chr5_147657804_147658307",]

#Maml3 +360kb 
peakGene[as.character(peakGene$geneName)=="Maml3",]
peakGene[as.character(peakGene$peakID)=="chr3_51744403_51744904",]


