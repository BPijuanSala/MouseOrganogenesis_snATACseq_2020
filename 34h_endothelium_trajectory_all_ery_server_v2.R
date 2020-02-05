###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Find dynamically accessible regions in snATAC-seq
###############################################################################


library(ggplot2)

#========functions=================
#These functions were developed for Pijuan-Sala et al ., Nature, 2019.

#This function determines the accessibility trajectories for clustering
get_predict = function(count, dpt,
                       model_mode = c("loess", "quadratic"),
                       predict_mode = c("cells", "even")){
  if(model_mode[1] == "loess"){
    model = loess(count ~ dpt)
    newdat = data.frame(dpt = dpt)
    if(predict_mode[1] == "even"){
      newdat = data.frame(dpt = seq(from = min(dpt), to = max(dpt), length.out = 1000))
    }
  } else if(model_mode[1] == "quadratic"){
    dpt2 = dpt^2
    model = lm(count ~ dpt + dpt2)
    newdat = data.frame(dpt = dpt, dpt2 = dpt2)
    if(predict_mode[1] == "even"){
      newdat = data.frame(dpt = seq(from = min(dpt), to = max(dpt), length.out = 1000),
                          dpt2 = seq(from = min(dpt), to = max(dpt), length.out = 1000)^2)
    }
  }
  pred = predict(model, newdata = newdat)#, dpt2 = coords^2))
  return(pred)
}
#This function identtifies genes that vary, and performs clustering
get_gene_clusters = function(sce, dpt, seed = 42, predict_even = TRUE, normalise = TRUE, model_mode = c("loess", "quadratic")){
  require(dynamicTreeCut)
  set.seed(seed)
  #find genes that vary along DPT
  dpt2 = dpt^2
  retain = sapply(1:nrow(sce), function(x){
    expr = sce[x,]
    mod = lm(expr ~ dpt + dpt2)
    base = lm(expr ~ 1)
    # dAIC = AIC(mod) - AIC(base)
    # return(dAIC < -25)
    anova = anova(base,mod)
    return(anova$`Pr(>F)`[2])
  })
  retain = p.adjust(retain, "fdr") < 0.05
  sce = sce[retain,]
  #fit the genes that vary
  predict_mode = ifelse(predict_even, "even", "cells")
  predicts = lapply(1:nrow(sce), function(x) get_predict(as.numeric(sce[x,]),
                                                         dpt,
                                                         predict_mode = predict_mode,
                                                         model_mode = model_mode[1]))
  predicts = do.call(cbind, predicts)
  if(normalise){
    #shift minimum value to 0
    predicts = sweep(predicts,
                     2,
                     matrixStats::colMins(predicts),
                     "-")
    #scale to max 1
    predicts = sweep(predicts,
                     2,
                     matrixStats::colMaxs(predicts),
                     "/")
  }
  colnames(predicts) = rownames(sce)
  #cluster
  dm = sqrt( (1- cor(predicts, method = "pearson")) /2)
  # dm = as.matrix(dist(t(predicts)))
  hclust = hclust(as.dist(dm), method = "average")
  clusts = cutreeDynamic(dendro = hclust, minClusterSize = 50, distM = dm)
  names(clusts) = colnames(predicts)
  return(list(clusters = clusts, hclust = hclust, predictions = predicts, dpt = dpt, expression = sce))
}


####
get_patterns = function(sce, dpt, region,seed = 42, predict_even = TRUE, normalise = TRUE, model_mode = c("loess", "quadratic")){
  set.seed(seed)
  #find genes that vary along DPT

  sce = sce[region,]
  #fit the genes that vary
  predict_mode = ifelse(predict_even, "even", "cells")
  if (length (region > 1)){
    
    predicts = lapply(1:nrow(sce), function(x) get_predict(as.numeric(sce[x,]),
                                                           dpt,
                                                           predict_mode = predict_mode,
                                                           model_mode = model_mode[1]))
    predicts = do.call(cbind, predicts)
    
    
    
    if(normalise){
      #shift minimum value to 0
      predicts = sweep(predicts,
                       2,
                       matrixStats::colMins(predicts),
                       "-")
      #scale to max 1
      predicts = sweep(predicts,
                       2,
                       matrixStats::colMaxs(predicts),
                       "/")
    }
    colnames(predicts) = rownames(sce)
    
  } else {
    predicts = get_predict(as.numeric(sce[region,]),
                dpt,
                predict_mode = predict_mode,
                model_mode = model_mode[1])
    
    
    if(normalise){
      #shift minimum value to 0
      predicts = predicts/max(predicts)
   
    }

    
    
  }
  
  return(predicts)
}
####
plot_cluster_averages = function(get_gene_clusters_out, predict_even = FALSE){
  input = get_gene_clusters_out
  means = sapply(1:length(unique(input$clusters)), function(x){
    rowMeans(input$predictions[,input$clusters==x])
  })
  pdf = melt(means)
  names(pdf) = c("cellnum", "cluster", "val")
  
  if(predict_even){
    pdf$dpt = pdf$cellnum
  } else{
    pdf$dpt = input$dpt[pdf$cellnum]
  }
  
  pdf$sd = sapply(1:nrow(pdf), function(x){
    sd(input$predictions[pdf$cellnum[x],
                         input$clusters == pdf$cluster[x]])
  })
  plots = lapply(1:length(unique(pdf$cluster)), function(x){
    ggplot(pdf[pdf$cluster == x,], aes(x= dpt, y = val)) +
      geom_ribbon(mapping = aes(ymin = val-sd, ymax = val + sd), fill = "grey80") +
      geom_line() + theme_bw() + 
      theme(#axis.text = element_blank(),  #axis.ticks = element_blank()
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      ggtitle(paste0(x, " (n=", sum(input$cluster==x),")" ))
  })
  return(plot_grid(plotlist = plots))
}


#========================
#Get dynamic patterns
#========================

wd = "/path/to/directory/sample_pooled_preprocess_revision1/"
library(moduleColor)
library(ggplot2)
bluePal <- c("#BFBFBF","#6495ED","#000000")

#Read metadata
meta = read.table(paste0(wd,"14_visualisation/data/snATACseq_embryo_revision1_metadata_afterAllQC_passQC_UMAP_cluster.txt"),
                  sep="\t",header=T)
meta$final_clusters = paste0("cluster_",meta$final_clusters)


meta_endo = meta[as.character(meta$final_clusters)%in%c("cluster_2","cluster_9","cluster_13"),]
rownames(meta_endo)=meta_endo$barcode

cellNames<-read.csv(paste0(wd,"18_endothelium_analysis/data/snATACseq_embryo_revision1_endothelium_eryAl_cellNames.csv"))
cellNames2<-read.csv(paste0(wd,"18_endothelium_analysis/data/snATACseq_embryo_revision1_endothelium_cellNames.csv"))

faScanpy<-read.csv(paste0(wd,"18_endothelium_analysis/data/snATACseq_embryo_revision1_endothelium_eryAl_FA.csv"))

clustScanpy<-read.csv(paste0(wd,"18_endothelium_analysis/data/snATACseq_embryo_revision1_endothelium_eryAl_louvain.csv"))



meta_endo[as.character(cellNames$X0),"FA_X"] = faScanpy$X0
meta_endo[as.character(cellNames$X0),"FA_Y"] = faScanpy$X1

meta_endo[as.character(cellNames$X0),"louvain_subclust"] = paste0("subclust_",clustScanpy$louvain)

#Read matrix
m <- Matrix::readMM(paste0(wd, '11_matrix_afterClusterQC/embryo_revision1_allPeaks_afterClusterPeak.mat.bin.mtx'))*1


cat(dim(m))

cat("\n reading barcodes\n")

barcodes = as.vector(read.csv(paste0(wd,"04_barcodeStats/embryo_revision1.xgi"),
                              sep="\n",header=F)[,1])   


cat(dim(barcodes))
colnames(m) <- barcodes
cat("\n reading bed file\n")

regions_frame = read.table(paste0(wd,"10_peaks_all/peaks_celltypes_extended500bp_qvalThreshold_overlapped_TSSGeneralPeaks_merged_noblacklist.bed"),
                           header=F, sep="\t")
colnames(regions_frame) <- c('Chr', 'Start', 'End')

rownames(m) <- paste0(regions_frame[,1], '_', regions_frame[,2], '_', regions_frame[,3])


#####======================================
# ALLANTOIS TO ENDOTHELIUM DPT
#####======================================
meta_endo1 = meta_endo[meta_endo$louvain_subclust%in%c("subclust_0","subclust_4","subclust_7"),]
rownames(meta_endo1)=meta_endo1$barcode
dpt<-read.csv(paste0(wd,"18_endothelium_analysis/data/snATACseq_embryo_revision1_endothelium_eryAl_dpt_sub_DPTtraj.csv"))
cellNames3<-read.csv(paste0(wd,"18_endothelium_analysis/data/snATACseq_embryo_revision1_endothelium_eryAl_cellNames_sub_DPTtraj.csv"))

rownames(dpt) = cellNames3$X0

meta_endo1[rownames(dpt) ,"dpt"] = dpt$dpt_pseudotime


barcodesPass = rownames(meta_endo1)

dat_endo = m[,barcodesPass]



dataPDTprob <- dat_endo[,order(meta_endo1$dpt,decreasing = F)]

#remove regions with where <10 cells present accessibility
dataPDTprob = dataPDTprob[Matrix::rowSums(dataPDTprob)>10,]
pdt = meta_endo1$dpt[order(meta_endo1$dpt,decreasing = F)]


clusters = get_gene_clusters(dataPDTprob, 
                                 pdt, 
                                 seed = 42, 
                                 predict_even = TRUE, 
                                 normalise = TRUE,
                                 model_mode = "loess")

library(reshape)
library(cowplot)
clusters_plot = plot_cluster_averages(clusters, predict_even = TRUE)
clusters_plot
ggsave(clusters_plot,
       file = paste0(wd,"18_endothelium_analysis/clusters_as_in_gutPijuanSala.pdf"),
       width = 8, height = 8)


clust.genes = clusters$clusters
for (i in 1:length(table(clust.genes))){
  c = names(table(clust.genes))[i]
  coord = names(clust.genes)[clust.genes == c]
  write.table(coord,file=paste0(wd,"18_endothelium_analysis/data/DPT_clusters_asGut/snATACseq_endothelium_eryAl_DPT_cluster_",c,".txt"),
              sep="\n",quote=F,row.names=F,col.names=F)
  
}

save(clusters,file=paste0(wd,"18_endothelium_analysis/data/DPT_clusters_asGut/snATACseq_endothelium_eryAl_DPT_clusters.rda"))
save(clusters_plot,file=paste0(wd,"18_endothelium_analysis/data/DPT_clusters_asGut/snATACseq_endothelium_eryAl_DPT_clusters_plot.rda"))


##==========================================
####LOCALLY
##==========================================

#========================
#Read data
#========================
clust_endo_colours = c(
  "0"="#f9a602",#Endothelium
  "1"="#8d021f",#Ery
  "2"="#933a16",
  #"2"="#b43757",#Ery
  "8"="#5e1914",#Ery
  "4"="#c06c84",#transition AL
  "5"="#260805",#Ery
  "6"="#ef820d",#HE
  "7"="#81007f",#AL
  
  "3"="#bf0a30",#Ery
  "9"="#7852a9",#AL
  "10"="#b200ed",#AL
  
  "11"="#ea3c53",#Ery
  "12"="#a45a52",#Ery
  "13"="#6f2da8",#AL
  "14"= "#fa8072",#Ery
  "15"="#ff0800",#transition blood
  "16"="#ed2939"
  
)

clust_endo = c(
  "0"="EC1",#Endothelium
  "1"="Ery1",#Ery
  "2"="Ery2",
  #"2"="#b43757",#Ery
  "8"="Ery3",#Ery
  "4"="Al_EC",#transition AL
  "5"="Ery4",#Ery
  "6"="EC2",#HE
  "7"="Al1",#AL
  
  "3"="Ery5",#Ery
  "9"="Al2",#AL
  "10"="Al3",#AL
  
  "11"="Ery6",#Ery
  "12"="Ery7",#Ery
  "13"="Al4",#AL
  "14"= "Ery8",#Ery
  "15"="EC_Haem",#transition blood
  "16"="Ery9"
  
)
library(ggplot2)
wd <- "/path/to/directory/"
load(file=paste0(wd,"PhD_BPS53/data/revision1/18_endothelium_analysis/cluster_DPT/snATACseq_endothelium_eryAl_DPT_clusters_plot.rda"))
load(file=paste0(wd,"PhD_BPS53/data/revision1/18_endothelium_analysis/cluster_DPT/snATACseq_endothelium_eryAl_DPT_clusters.rda"))
library(reshape2)
library(cowplot)
clusters_plot = plot_cluster_averages(clusters, predict_even = TRUE)
clusters_plot

#ggsave(clusters_plot,
#       file = paste0(wd,"PhD_BPS53/plots/revision1/18_endothelium_analysis/clusters_as_in_gutPijuanSala_endoToAll.pdf"),
#       width = 9, height = 7)


# Read metadata

meta = read.table(paste0(wd,"PhD_BPS53/data/revision1/Supplementary_Table5_metadata_passQC_topics.txt"),
                  sep="\t",header=T)
meta$final_clusters = paste0("cluster_",meta$final_clusters)

meta_endo = meta[as.character(meta$final_clusters)%in%c("cluster_2","cluster_9","cluster_13"),]
rownames(meta_endo)=meta_endo$barcode

cellNames<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/18_endothelium_analysis/snATACseq_embryo_revision1_endothelium_eryAl_cellNames.csv"))
cellNames2<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/18_endothelium_analysis/snATACseq_embryo_revision1_endothelium_cellNames.csv"))

faScanpy<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/18_endothelium_analysis/snATACseq_embryo_revision1_endothelium_eryAl_FA.csv"))

clustScanpy<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/18_endothelium_analysis/snATACseq_embryo_revision1_endothelium_eryAl_louvain.csv"))



meta_endo[as.character(cellNames$X0),"FA_X"] = faScanpy$X0
meta_endo[as.character(cellNames$X0),"FA_Y"] = faScanpy$X1

meta_endo[as.character(cellNames$X0),"louvain_subclust"] = paste0("subclust_",clustScanpy$louvain)


#Subset clusters for DPT
meta_endo1 = meta_endo[meta_endo$louvain_subclust%in%c("subclust_0","subclust_4","subclust_7"),]
rownames(meta_endo1)=meta_endo1$barcode
dpt<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/18_endothelium_analysis/snATACseq_embryo_revision1_endothelium_eryAl_dpt_sub_DPTtraj.csv"))
cellNames3<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/18_endothelium_analysis/snATACseq_embryo_revision1_endothelium_eryAl_cellNames_sub_DPTtraj.csv"))

rownames(dpt) = cellNames3$X0

meta_endo1[rownames(dpt) ,"dpt"] = dpt$dpt_pseudotime

meta_endo1[rownames(dpt) ,"dpt_x1000"] = dpt$dpt_pseudotime*1000

# Plot FDG
tiff(paste0(wd,"PhD_BPS53/plots/revision1/18_endothelium_analysis/FDG_endothelium_ery_al.tiff"),
     width=7,height=7,res=300,units="in")
plot(meta_endo$FA_X,meta_endo$FA_Y,pch=20,axes=F,ylab="",xlab="",
     #col=clust_endo_colours[17])
     col=unname(clust_endo_colours[as.character(gsub("subclust_","",meta_endo$louvain_subclust))]))
dev.off()

# Plot UMAP
tiff(paste0(wd,"PhD_BPS53/plots/revision1/18_endothelium_analysis/UMAP_endothelium_ery_al.tiff"),
     width=7,height=7,res=300,units="in")
plot(meta_endo$umap_X,meta_endo$umap_Y,pch=20,axes=F,ylab="",xlab="",
     #col=clust_endo_colours[17])
     col=unname(clust_endo_colours[as.character(gsub("subclust_","",meta_endo$louvain_subclust))]))
dev.off()

tiff(paste0(wd,"PhD_BPS53/plots/revision1/18_endothelium_analysis/UMAP_endothelium_ery_al_clusters.tiff"),
     width=7,height=7,res=300,units="in")
plot(meta_endo$umap_X,meta_endo$umap_Y,pch=20,axes=F,ylab="",xlab="",
     #col=clust_endo_colours[17])
     col=unname(all_colours[as.character(meta_endo$ann)]))
dev.off()




names(clust_endo_colours) = clust_endo
all_colours_sub =all_colours[as.character(unique(meta_endo$ann))]
p = ggplot(meta_endo, aes(x=umap_X, y=umap_Y))
plot = p + geom_point(aes(col = factor(meta_endo$ann))) + 
  scale_color_manual(values = unname(all_colours_sub[order(names(all_colours_sub))]), 
                     labels = names(all_colours_sub)[order(names(all_colours_sub))], 
                     drop = FALSE, name = "")  +
  theme_void() +
  theme(legend.text=element_text(size=rel(1.1))) +
  guides(colour = guide_legend(override.aes = list(size=7))) + ggtitle("")
print(plot)
ggsave(filename=paste0(wd,"PhD_BPS53/plots/revision1/UMAP_AL_EC_HAEM_celltypes_website.png"),
       plot=plot,width = 7, height = 6,
       device="png")



p = ggplot(meta_endo, aes(x=umap_X, y=umap_Y))
plot = p + geom_point(aes(col = factor(meta_endo$al_haem_endo_clusters))) + 
  scale_color_manual(values = unname(clust_endo_colours[order(names(clust_endo_colours))]), 
                     labels = names(clust_endo_colours)[order(names(clust_endo_colours))], 
                     drop = FALSE, name = "")  +
  theme_void() +
  theme(legend.text=element_text(size=rel(1.1))) +
  guides(colour = guide_legend(override.aes = list(size=7))) + ggtitle("")
print(plot)
ggsave(filename=paste0(wd,"PhD_BPS53/plots/revision1/UMAP_AL_EC_HAEM_website.png"),
       plot=plot,width = 7, height = 6,
       device="png")





# Plot gradient of pseudotime
meta_endo$DPT_together = rep(0,nrow(meta_endo))
meta_endo[rownames(meta_endo1),"DPT_together"] = meta_endo1$dpt

tiff(paste0(wd,"PhD_BPS53/plots/revision1/18_endothelium_analysis/UMAP_endothelium_ATACseq_meso_al_DPT_together.tiff"),
     width=7,height=7,res=300,units="in")
p = ggplot(meta_endo, aes(x=umap_X, y=umap_Y))
plot = p + geom_point(aes(col = DPT_together))+
  theme_void() +
  
  scale_fill_gradient2(low = "orange", mid = "#BFBFBF",
                       high = "#0664cf", midpoint = 0, space = "Lab",
                       na.value = "white", 
                       guide = "colourbar", aesthetics = "colour") +
  theme(legend.title = element_blank(),legend.position = "none") + 
  ggtitle("")
print(plot)
dev.off()

tiff(paste0(wd,"PhD_BPS53/plots/revision1/18_endothelium_analysis/FDG_endothelium_ATACseq_meso_al_DPT_together_withLegend.tiff"),
     width=7,height=7,res=300,units="in")
p = ggplot(meta_endo, aes(x=FA_X, y=FA_Y))
plot = p + geom_point(aes(col = DPT_together))+
  theme_void() +
  
  scale_fill_gradient2(low = "orange", mid = "#BFBFBF",
                       high = "#0664cf", midpoint = 0, space = "Lab",
                       na.value = "white", 
                       guide = "colourbar", aesthetics = "colour") +
  #theme(legend.title = element_blank(),legend.position = "none") + 
  ggtitle("")
print(plot)
dev.off()


#Plot cells along pseudotime with cluster colours
pdf(paste0(wd,"PhD_BPS53/plots/revision1/18_endothelium_analysis/DPT_allantois_subclusterscelltype.pdf"),
   width=5,height=2.7)
plot(x= meta_endo1[ ,"dpt_x1000"],
     y=as.numeric(as.factor(meta_endo1$louvain_subclust)),
     col=alpha(clust_endo_colours[gsub("subclust_","",meta_endo1$louvain_subclust)], 0.4), 
     ylab="",xlab="",pch=16,cex=1.2,ylim=c(0,3.2)
     )
dev.off()


##Patterns bound by Etv2
clust_num = table(clusters$clusters)
clust_etv2 = c(162,207,201,0,0,103,150,17,0,19,58,1)
etv2_ratio = sort(clust_etv2/clust_num)

pdf(paste0(wd,"PhD_BPS53/plots/revision1/18_endothelium_analysis/Patterns_bound_by_etv2.pdf"),
   width=5,height=4)

barplot(etv2_ratio,col="#0869a1",ylim=c(0,1))
dev.off()