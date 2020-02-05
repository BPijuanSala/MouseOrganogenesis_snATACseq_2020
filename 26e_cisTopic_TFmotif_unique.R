###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  October 2019
##  DESCRIPTION: Script to define regions uniquely contributing to each topic and their
##                TF enrichment.
###############################################################################


###===============================
## On the SERVER
###===============================
library(plyr)
suppressWarnings(library(cisTopic))

library(GenomicRanges)
library(rtracklayer)
library(Matrix)


wd = "/path/to/directory/sample_pooled_preprocess_revision1/"


####===============================================
##    Run Rcistarget as usual pipeline
####===============================================
load(file=paste0(wd,"13_cisTopic/data/cisTopicObject_50_100_afterReadsInPeaksQC_24.rda"))



clustScanpy<-read.csv(paste0(wd,"14_visualisation/data/snATACseq_embryo_revision1_reads_in_peaks_above_24_Louvain.csv"))
rownames(clustScanpy)<-clustScanpy[,1]
densityClust <- paste0("subclust_",clustScanpy$louvain)
densityClust <- as.data.frame(densityClust)
rownames(densityClust) <-rownames(clustScanpy)
colnames(densityClust) <- 'densityClust'
densityClust[,1] <- as.factor(densityClust[,1])

cisTopicObject <- addCellMetadata(cisTopicObject, densityClust)



# To use RcisTarget, we need to liftover the coordinates to the mm9 assemble
library(R.utils)
cat("\n downloading liftover\n")
url <- "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm9.over.chain.gz"
mm10Tomm9.chain <- "mm10Tomm9.over.chain"
download.file(url, destfile = paste0(mm10Tomm9.chain, ".gz"))
gunzip(paste0(mm10Tomm9.chain, ".gz"))

# Import chain file
mm10Tomm9.chain  <- import.chain(mm10Tomm9.chain)

# Obtain liftOver dictionary (as list)
mm10_coord <- cisTopicObject@region.ranges
cat("\n liftover\n")

mm10_to_mm9_list <- liftOver(mm10_coord, mm10Tomm9.chain)



# Run GREAT based on liftover to mm9 coordinates
cat("\n convert to binarised\n")
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)

cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)

cisTopicObject <- binarizedcisTopicsToCtx(cisTopicObject, liftOver=mm10_to_mm9_list, genome='mm9')

cat("\n Rscistarget\n")


pathToFeather <- "/path/to/directory/bioinformatic_resources/databases_motif/mm9-regions-9species.all_regions.mc9nr.feather"
cisTopicObject <- topicsRcisTarget(cisTopicObject, genome='mm9', 
                                   pathToFeather, 
                                   reduced_database=FALSE, 
                                   nesThreshold=3, rocthr=0.005, 
                                   maxRank=20000, nCores=10)


cat("\n save\n")

save(cisTopicObject,file=paste0(wd,"13_cisTopic/data/cisTopicObject_50_100_afterReadsInPeaksQC_24_Rcistarget.rda"))


####=================================================================================
##    Run Rcistarget on the regions uniquely contributing to each topic
####=================================================================================
load(file=paste0(wd,"13_cisTopic/data/cisTopicObject_50_100_afterReadsInPeaksQC_24_Rcistarget.rda"))

##Get only regions uniquely contributing to each topic
topicsList2 <- cisTopicObject@binarized.regions.to.Rct
topics = unique(unlist(topicsList2))

top_counts = data.frame(
  coord = topics,
  topic = rep(NA,length(topics))
)


topicsList= topicsList2
for (i in 1:100){
  tab = topicsList[[i]]
  top_counts$topic[match(tab,top_counts$coord)] = paste0(top_counts[match(tab,top_counts$coord),"topic"],"","Topic",i,";")
}
top_counts$topic = gsub("NAT","T",top_counts$topic)
top_counts$topic_stringent = top_counts$topic
top_counts$topic_stringent[grep(";Topic",top_counts$topic_stringent)] = "Nonspecific"

for (i in 1:100){
  topicsList2[[i]] = as.character(top_counts$coord)[as.character(top_counts$topic_stringent)==paste0("Topic",i,";")]
}
topicsList2[[100]]
names(topicsList2)
for (i in 1:100){
  if(length(topicsList2[[i]])==0){topicsList2[[i]] <- NULL}
}

# Run Rcistarget
pathToFeather <- "/path/to/directory/bioinformatic_resources/databases_motif/mm9-regions-9species.all_regions.mc9nr.feather"

cisTopicObject@binarized.regions.to.Rct = topicsList2
cisTopicObject <- topicsRcisTarget(cisTopicObject, genome='mm9', 
                                   pathToFeather, 
                                   reduced_database=FALSE, 
                                   nesThreshold=3, rocthr=0.005, 
                                   maxRank=20000, nCores=1)

RcistargetOutput = cisTopicObject@binarized.RcisTarget
names(RcistargetOutput) = names(topicsList2)
for (i in 1:length(RcistargetOutput)){
  dat = RcistargetOutput[[i]]
  nameT = names(RcistargetOutput)[i]
  write.table(dat,file=paste0(wd,"13_cisTopic/data/cisTopicObject_50_100_afterReadsInPeaksQC_24_Rcistarget_unique_",nameT,".txt"),
              sep="\t",quote=F,row.names=T,col.names=T)
}
save(RcistargetOutput,file=paste0(wd,"13_cisTopic/data/cisTopicObject_50_100_afterReadsInPeaksQC_24_Rcistarget_unique_object_TF.rda"))

##--------------------------------------
#Computed LOCALLY
##--------------------------------------
library(ggplot2)
wd <- "/path/to/directory/"

#########TF enrichment analysis on Topics
load(paste0(wd,"PhD_BPS53/data/revision1/13_cisTopic/topics_Rcistarget/cisTopicObject_50_100_afterReadsInPeaksQC_24_Rcistarget_unique_object_TF.rda"))


#NMPs
Topic38 <- RcistargetOutput[["Topic38"]]
Topic51 <- RcistargetOutput[["Topic51"]]
Topic100 <- RcistargetOutput[["Topic100"]]


DT::datatable(Topic51[,(colnames(Topic51)%in% c("enrichedRegions", "TF_lowConf"))==FALSE], 
              escape = FALSE, filter="top", options=list(pageLength=100))
motifs = DT::datatable(Topic100[,(colnames(Topic38)%in% c("enrichedRegions", "TF_lowConf"))==FALSE], 
              escape = FALSE, filter="top", options=list(pageLength=100))
motifs2 = RcisTarget::addLogo(motifs)

d = merge(merge(Topic38,Topic100,by="motif",all=T),Topic51,by="motif",all=T)

e = d[,c('motif','NES.x',#'AUC.x','TF_highConf.x',
         'NES.y',#'AUC.y',#'TF_highConf.y',
         'NES',#'AUC',
         'TF_highConf'
         )]

#After checking the results, define the motifs to plot
motifs = c(
'cisbp__M5312',
  
  'taipale_cyt_meth__CDX2_NYAATAAAN_eDBD',
  "taipale_cyt_meth__HOXC11_RGYAATAAAAN_FL_meth",
  'taipale_cyt_meth__HOXA11_RGYAATAAAAN_eDBD_meth',
  'hocomoco__EVX2_MOUSE.H11MO.0.A',
  'transfac_pro__M08973',
  'idmmpmm__D')


d[motifs,c("TF_highConf","logo")]
f = e[order(e$NES,decreasing=T),]
f[1:20,]

sel = f[match(motifs,f$motif),]
sel



#Plot motif scores in heatmap
heatmapRedYelBlue <- c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090",
                       "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")

palette <- colorRampPalette(rev(heatmapRedYelBlue))

rownames(sel) = gsub(" .*","",sel$TF_highConf)
sel2 = as.matrix(sel[,c("NES.x","NES.y","NES")])
sel2[is.na(sel2)] = 0


pdf(paste0(wd,"PhD_BPS53/plots/revision1/heatmap_NMP_topics_TFs.pdf"),
    width=4,height=7)
gplots::heatmap.2(sel2, trace="none", 
                  col=palette, 
                  Colv = F,
                  Rowv = T, 
                  #ColSideColors = c, 
                  #RowSideColors = cols, 
                  dendrogram = "none", density.info = 'none',
                  #labRow = g,
                  key=T#labCol=""
)
dev.off()


