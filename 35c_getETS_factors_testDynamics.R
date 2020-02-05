###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  October 2019
##  DESCRIPTION: Test dynamics of ETS factors
###############################################################################

library(ggplot2)

#========functions=================
#These functions were developed for Pijuan-Sala et al ., Nature, 2019.

#This function determines the expression trajectories for clustering
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
  clusts = cutreeDynamic(dendro = hclust, minClusterSize = 2, distM = dm)
  names(clusts) = colnames(predicts)
  return(list(clusters = clusts, hclust = hclust, predictions = predicts, dpt = dpt, expression = sce))
}
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
      geom_line() +
      theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
      ggtitle(paste0(x, " (n=", sum(input$cluster==x),")" ))
  })
  return(plot_grid(plotlist = plots))
}
#========================

ETS = anSeq::getGeneID(c("Elf1","Elf2","Elf4",
                   "Gabpa",
                   "Erg","Fli1","Fev",
                   "Erf","Etv3",
                   "Elf3","Elf5","Ese3",
                   "Ets1","Ets2",
                   "Spdef",
                   "Etv4","Etv5","Etv1",
                   "Etv2",
                   "Spi1","Elk4","Elk3",
                   "Etv6","Etv7"
                   ))[[2]]

anSeq::getGeneID(c("Elf1","Elf2","Elf4",
                   "Gabpa",
                   "Erg","Fli1","Fev",
                   "Erf","Etv3",
                   "Elf3","Elf5","Ese3",
                   "Ets1","Ets2",
                   "Spdef",
                   "Etv4","Etv5","Etv1",
                   "Etv2",
                   "Spi1","Elk4","Elk3",
                   "Etv6","Etv7"
))[[1]]

ETS = c("ENSMUSG00000036461", "ENSMUSG00000037174", "ENSMUSG00000031103" ,"ENSMUSG00000008976",
               "ENSMUSG00000040732" ,"ENSMUSG00000016087", "ENSMUSG00000055197", "ENSMUSG00000040857",
               "ENSMUSG00000003382", "ENSMUSG00000003051", "ENSMUSG00000027186", "ENSMUSG00000032035",
               "ENSMUSG00000022895" ,"ENSMUSG00000024215" ,"ENSMUSG00000017724" ,"ENSMUSG00000013089",
               "ENSMUSG00000004151" ,"ENSMUSG00000006311", "ENSMUSG00000002111" ,"ENSMUSG00000026436",
               "ENSMUSG00000008398" ,"ENSMUSG00000030199",
               
               'ENSMUSG00000031871','ENSMUSG00000020717')
ETS_names = c("Elf1",  "Elf2" , "Elf4" , "Gabpa", "Erg" ,  "Fli1" , "Fev"  , "Erf" ,  "Etv3"  ,"Elf3", 
             "Elf5" , "Ets1" , "Ets2" , "Spdef", "Etv4" , "Etv5"  ,"Etv1" , "Etv2" , "Spi1" , "Elk4",
             "Elk3" , "Etv6","Cdh5",'Pecam1')

wd = "/path/to/directory/sample_pooled_preprocess_revision1/"
library(moduleColor)
library(ggplot2)
bluePal <- c("#BFBFBF","#6495ED","#000000")

#Read metadata RNAseq
metadataFilename ='/path/to/directory/scRNAseq_PijuanSala_2019/data/PijuanSalaEtAl_SupplementaryTable3_revision.txt'

meta = read.csv(metadataFilename,
                  sep="\t",header=T)

rownames(meta)=meta$cell

cellNames<-read.csv(paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_embryo_revision1_endothelium_scRNAseq_cellNames.csv"))
meta_endo = meta[cellNames$X0,]
faScanpy<-read.csv(paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_embryo_revision1_endothelium_scRNAseq_FA.csv"))

meta_endo[as.character(cellNames$X0),"FA_X"] = faScanpy$X0
meta_endo[as.character(cellNames$X0),"FA_Y"] = faScanpy$X1

clustScanpy<-read.csv(paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_embryo_revision1_endothelium_scRNAseq_Louvain.csv"))

meta_endo[as.character(cellNames$X0),"louvain_subclust"] = paste0("subclust_",clustScanpy$louvain)


#Allantois trajectory
cellNames_al<-read.csv(paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_embryo_revision1_endothelium_scRNAseq_allantois_DPTtraj_cellNames.csv"))
dpt_al<-read.csv(paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_embryo_revision1_endothelium_scRNAseq_allantois_DPTtraj_dpt.csv"))

meta_endo$DPT_al = rep(NA,nrow(meta_endo))
meta_endo[as.character(cellNames_al$X0),"DPT_al"] = dpt_al$dpt_pseudotime

#Mixed mesoderm trajectory
cellNames_mixed<-read.csv(paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_embryo_revision1_endothelium_scRNAseq_MixedMeso_DPTtraj_cellNames.csv"))
dpt_mixed<-read.csv(paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_embryo_revision1_endothelium_scRNAseq_MixedMeso_DPTtraj_dpt.csv"))

meta_endo$DPT_meso = rep(NA,nrow(meta_endo))
meta_endo[as.character(cellNames_mixed$X0),"DPT_meso"] = dpt_mixed$dpt_pseudotime

#Read counts matrix
m<- Matrix::readMM(paste0(wd, '21_endothelium_scRNAseq/data/countsData_endothelium_logged.mtx'))*1


cat(dim(m))

cat("\n reading barcodes\n")

cells = as.vector(read.csv(paste0(wd,"21_endothelium_scRNAseq/data/cells_endothelium.txt"),
                              sep="\n",header=F)[,1])   
genes = as.vector(read.csv(paste0(wd,"21_endothelium_scRNAseq/data/genes_endothelium.txt"),
                              sep="\n",header=F)[,1])   

colnames(m) <- cells
rownames(m) <- genes

##=========================================================================
## Get dynamically expressed ETS factors along the allantoic trajectory
##=========================================================================
dat_al = m[,rownames(meta_endo)[is.na(meta_endo$DPT_al)==F]]
meta_al = meta_endo[rownames(meta_endo)[is.na(meta_endo$DPT_al)==F],]

dataPDTprob <- dat_al[,order(meta_al$DPT_al,decreasing = F)]
dataPDTprob = dataPDTprob[rownames(dataPDTprob)%in%ETS,]
dataPDTprob = dataPDTprob[matrixStats::rowVars(as.matrix(dataPDTprob))>0.15,]#Filter for highly variable


pdt = meta_al$DPT_al[order(meta_al$DPT_al,decreasing = F)]

#normalised
clusters = get_gene_clusters(dataPDTprob, 
                             pdt, 
                             seed = 42, 
                             predict_even = TRUE, 
                             normalise = TRUE,
                             model_mode = "loess")
save(clusters,file=paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_endothelium_scRNAseq_ETS_allantois_dpt_varAbovepoint15_clusters_Normalised.rda"))

#not normalised
clusters_notnorm = get_gene_clusters(dataPDTprob, 
                             pdt, 
                             seed = 42, 
                             predict_even = TRUE, 
                             normalise = F,
                             model_mode = "loess")
save(clusters_notnorm,file=paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_endothelium_scRNAseq_ETS_allantois_varAbovepoint15_dpt_clusters_notNormalised.rda"))

# ..cutHeight not given, setting it to 0.708  ===>  99% of the (truncated) height range in dendro.
#..done.


library(reshape)
library(cowplot)
clusters_plot = plot_cluster_averages(clusters, predict_even = TRUE)
clusters_plot
ggsave(clusters_plot,
       file = paste0(wd,"21_endothelium_scRNAseq/plots/snATACseq_endothelium_scRNAseq_ETS_allantois_dpt_clusters_as_in_gutPijuanSala_varAbovepoint15.pdf"),
       width = 8, height = 8)


#save(clusters,file=paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_endothelium_scRNAseq_ETS_allantois_dpt_clusters.rda"))
save(clusters_plot,file=paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_endothelium_scRNAseq_ETS_allantois_dpt_clusters_plot_varAbovepoint15.rda"))


##=========================================================================
## Get dynamically expressed ETS factors along the mixed mesoderm trajectory
##=========================================================================

dat_meso = m[,rownames(meta_endo)[is.na(meta_endo$DPT_meso)==F]]
meta_meso = meta_endo[rownames(meta_endo)[is.na(meta_endo$DPT_meso)==F],]

dataPDTprob <- dat_meso[,order(meta_meso$DPT_meso,decreasing = F)]
dataPDTprob = dataPDTprob[rownames(dataPDTprob)%in%ETS,]
dataPDTprob = dataPDTprob[matrixStats::rowVars(as.matrix(dataPDTprob))>0.15,]#Filter for highly variable

pdt = meta_meso$DPT_meso[order(meta_meso$DPT_meso,decreasing = F)]

#normalised
clusters = get_gene_clusters(dataPDTprob, 
                             pdt, 
                             seed = 42, 
                             predict_even = TRUE, 
                             normalise = TRUE,
                             model_mode = "loess")

save(clusters,file=paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_endothelium_scRNAseq_ETS_mesoderm_dpt_clusters_Normalised_varAbovepoint15.rda"))

#not normalised
clusters_notnorm = get_gene_clusters(dataPDTprob, 
                                     pdt, 
                                     seed = 42, 
                                     predict_even = TRUE, 
                                     normalise = F,
                                     model_mode = "loess")
save(clusters_notnorm,file=paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_endothelium_scRNAseq_ETS_mesoderm_dpt_clusters_notNormalised_varAbovepoint15.rda"))

library(reshape)
library(cowplot)
clusters_plot = plot_cluster_averages(clusters, predict_even = TRUE)
clusters_plot
ggsave(clusters_plot,
       file = paste0(wd,"21_endothelium_scRNAseq/plots/snATACseq_endothelium_scRNAseq_ETS_mesoderm_dpt_clusters_as_in_gutPijuanSala_varAbovepoint15.pdf"),
       width = 8, height = 8)


#save(clusters,file=paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_endothelium_scRNAseq_ETS_mesoderm_dpt_clusters.rda"))
save(clusters_plot,file=paste0(wd,"21_endothelium_scRNAseq/data/snATACseq_endothelium_scRNAseq_ETS_mesoderm_dpt_clusters_plot_varAbovepoint15.rda"))



#==============================
####LOCALLY
#==============================

all_colours = c(
  "Mixed mesoderm" =  "#DFCDE4",  
  
  "Allantois" = "#532C8A",
  "Blood progenitors 1" = "#f9decf",
  "Blood progenitors 2" = "#c9a997",

  "Endothelium" =  "#ff891c",                               
  "Erythroid1" =  "#C72228",                               
  "Erythroid2" =  "#EF4E22",
  "Erythroid3" =  "#f77b59",
                                   
  "Haematoendothelial progenitors" =  "#FBBE92"         
  
)

spectralPal = c(
  'E6.5'="#D53E4F",
  'E6.75'="#F46D43",
  'E7.0'="#FDAE61",
  'E7.5'="#FFFFBF",
  'E7.25'="#FEE08B",
  'E7.75'="#E6F598",
  'E8.0'="#ABDDA4",
  'E8.5'="#3288BD",
  'E8.25'="#66C2A5",
  'mixed_gastrulation'= "#A9A9A9"  
  
)
library(ggplot2)
wd <- "/path/to/directory/"

ETS = c("ENSMUSG00000036461", "ENSMUSG00000037174", "ENSMUSG00000031103" ,"ENSMUSG00000008976",
        "ENSMUSG00000040732" ,"ENSMUSG00000016087", "ENSMUSG00000055197", "ENSMUSG00000040857",
        "ENSMUSG00000003382", "ENSMUSG00000003051", "ENSMUSG00000027186", "ENSMUSG00000032035",
        "ENSMUSG00000022895" ,"ENSMUSG00000024215" ,"ENSMUSG00000017724" ,"ENSMUSG00000013089",
        "ENSMUSG00000004151" ,"ENSMUSG00000006311", "ENSMUSG00000002111" ,"ENSMUSG00000026436",
        "ENSMUSG00000008398" ,"ENSMUSG00000030199",
        
        'ENSMUSG00000031871','ENSMUSG00000020717')
ETS_names = c("Elf1",  "Elf2" , "Elf4" , "Gabpa", "Erg" ,  "Fli1" , "Fev"  , "Erf" ,  "Etv3"  ,"Elf3", 
              "Elf5" , "Ets1" , "Ets2" , "Spdef", "Etv4" , "Etv5"  ,"Etv1" , "Etv2" , "Spi1" , "Elk4",
              "Elk3" , "Etv6","Cdh5",'Pecam1')


#read metadata scRNAseq
metadataFilename =paste0(wd,'PhD_BPS32/release6/paper_resources/important_data//PijuanSalaEtAl_SupplementaryTable3_revision.txt')

meta = read.csv(metadataFilename,
                sep="\t",header=T)

rownames(meta)=meta$cell

cellNames<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/21_endothelium_scRNAseq/snATACseq_embryo_revision1_endothelium_scRNAseq_cellNames.csv"))
meta_endo = meta[as.character(cellNames$X0),]
faScanpy<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/21_endothelium_scRNAseq/snATACseq_embryo_revision1_endothelium_scRNAseq_FA.csv"))

meta_endo[as.character(cellNames$X0),"FA_X"] = faScanpy$X0
meta_endo[as.character(cellNames$X0),"FA_Y"] = faScanpy$X1

clustScanpy<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/21_endothelium_scRNAseq//snATACseq_embryo_revision1_endothelium_scRNAseq_Louvain.csv"))

meta_endo[as.character(cellNames$X0),"louvain_subclust"] = paste0("subclust_",clustScanpy$louvain)


#Read DPT
cellNames_al<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/21_endothelium_scRNAseq//snATACseq_embryo_revision1_endothelium_scRNAseq_allantois_DPTtraj_cellNames.csv"))
dpt_al<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/21_endothelium_scRNAseq//snATACseq_embryo_revision1_endothelium_scRNAseq_allantois_DPTtraj_dpt.csv"))

meta_endo$DPT_al = rep(NA,nrow(meta_endo))
meta_endo[as.character(cellNames_al$X0),"DPT_al"] = dpt_al$dpt_pseudotime

cellNames_mixed<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/21_endothelium_scRNAseq//snATACseq_embryo_revision1_endothelium_scRNAseq_MixedMeso_DPTtraj_cellNames.csv"))
dpt_mixed<-read.csv(paste0(wd,"PhD_BPS53/data/revision1/21_endothelium_scRNAseq//snATACseq_embryo_revision1_endothelium_scRNAseq_MixedMeso_DPTtraj_dpt.csv"))

meta_endo$DPT_meso = rep(NA,nrow(meta_endo))
meta_endo[as.character(cellNames_mixed$X0),"DPT_meso"] = dpt_mixed$dpt_pseudotime

meta_endo[ ,"DPT_al_x1000"] = meta_endo$DPT_al*1000
meta_endo[ ,"DPT_meso_x1000"] = meta_endo$DPT_meso*1000

clust_endo_RNA = c(
  "0"="#4f97a3",#Haem from mix
  "1"="#6f2da8",#AL
  "2"="#95c8d8",#Mixed
  "3"="#f9a602",#Endo
  
  "4"="#c06c84",#al transition
  "5"="#b200ed",#al
  "6"="#7ef9ff",#mixed
  "7"="#fa8072",#Endo-haem
  "8"="#3fe0d0",#Mix transition
  
  "9"="#588bae",#mixed
  "10"="#ef820d",#Endo
  
  "11"="#ea3c53"#haem
  
)

# Plot FA coloured by subclusters
tiff(paste0(wd,"PhD_BPS53/plots/revision1/21_endothelium_scRNAseq/FDG_endothelium_scRNAseq_meso_al.tiff"),
     width=7,height=7,res=300,units="in")
plot(meta_endo$FA_X,meta_endo$FA_Y,pch=20,axes=F,ylab="",xlab="",
     #col=clust_endo_colours[17])
     col=unname(clust_endo_RNA[as.character(gsub("subclust_","",meta_endo$louvain_subclust))]))
dev.off()

# Plot FA coloured by cell type
tiff(paste0(wd,"PhD_BPS53/plots/revision1/21_endothelium_scRNAseq/FDG_endothelium_scRNAseq_meso_al_celltypes_original.tiff"),
     width=7,height=7,res=300,units="in")
plot(meta_endo$FA_X,meta_endo$FA_Y,pch=20,axes=F,ylab="",xlab="",
     #col=clust_endo_colours[17])
     col=unname(all_colours[as.character(meta_endo$celltype)]))
dev.off()

# Plot FA coloured by stage
tiff(paste0(wd,"PhD_BPS53/plots/revision1/21_endothelium_scRNAseq/FDG_endothelium_scRNAseq_meso_al_celltypes_original_stage.tiff"),
     width=7,height=7,res=300,units="in")
plot(meta_endo$FA_X,meta_endo$FA_Y,pch=20,axes=F,ylab="",xlab="",
     #col=clust_endo_colours[17])
     col=unname(spectralPal[as.character(meta_endo$stage)]))
dev.off()

for (i in sort(unique(as.character(meta_endo$stage)))){
  plot(meta_endo$FA_X,meta_endo$FA_Y,pch=20,axes=F,ylab="",xlab="",main=i,
       #col=clust_endo_colours[17])
       col="gray")
  points(meta_endo$FA_X[as.character(meta_endo$stage)==i],
         meta_endo$FA_Y[as.character(meta_endo$stage)==i],
         pch=20,axes=F,ylab="",xlab="",
       #col=clust_endo_colours[17])
       col=unname(spectralPal[as.character(i)]))
}


##PSEUDOTIME PLOTS===================================================

#Combine DPT scores for plot
meta_endo$DPT_together = rep(0,nrow(meta_endo))
meta_endo$DPT_together[is.na(meta_endo$DPT_al)==F] = meta_endo$DPT_al[is.na(meta_endo$DPT_al)==F]
meta_endo$DPT_together[is.na(meta_endo$DPT_meso)==F] = -meta_endo$DPT_meso[is.na(meta_endo$DPT_meso)==F]
tiff(paste0(wd,"PhD_BPS53/plots/revision1/21_endothelium_scRNAseq/FDG_endothelium_scRNAseq_meso_al_DPT_together.tiff"),
     width=7,height=7,res=300,units="in")
p = ggplot(meta_endo, aes(x=FA_X, y=FA_Y))
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



tiff(paste0(wd,"PhD_BPS53/plots/revision1/21_endothelium_scRNAseq/FDG_endothelium_scRNAseq_meso_al_DPT_together_withLegend.tiff"),
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


#===================================================


#PLOTS CELLTYPES (dots) ALONG TRANSITION===================================================


pdf(paste0(wd,"PhD_BPS53/plots/revision1/21_endothelium_scRNAseq/scRNAseq_DPT_AL_to_endo_celltypes_points.pdf"),
    width=5,height=2.7)
tmp = meta_endo[is.na(meta_endo$DPT_al)==F,]
plot(x= tmp[ ,"DPT_al_x1000"],
     y=as.numeric(as.factor(tmp$louvain_subclust)),
     col=alpha(clust_endo_RNA[gsub("subclust_","",tmp$louvain_subclust)], 0.4), 
     ylab="",xlab="",pch=16,cex=1.2,ylim=c(0,3.2)
)
dev.off()


pdf(paste0(wd,"PhD_BPS53/plots/revision1/21_endothelium_scRNAseq/scRNAseq_DPT_AL_to_endo_celltypes_points_original.pdf"),
    width=5,height=2.7)
plot(x= tmp[ ,"DPT_al_x1000"],
     y=as.numeric(as.factor(tmp$louvain_subclust)),
     col=alpha(all_colours[as.character(tmp$celltype)], 0.4), 
     ylab="",xlab="",pch=16,cex=1.2,ylim=c(0,3.2)
)
dev.off()



pdf(paste0(wd,"PhD_BPS53/plots/revision1/21_endothelium_scRNAseq/scRNAseq_DPT_meso_to_endo_celltypes_points.pdf"),
    width=5,height=2.7)
tmp = meta_endo[is.na(meta_endo$DPT_meso)==F,]
tmp = tmp[order(tmp$DPT_meso_x1000,decreasing=F),]
plot(x= tmp[ ,"DPT_meso_x1000"],
     y=as.numeric(as.factor(tmp$louvain_subclust)),
     col=alpha(clust_endo_RNA[gsub("subclust_","",tmp$louvain_subclust)], 0.4), 
     ylab="",xlab="",pch=16,cex=1.2,ylim=c(0,5.2)
)
dev.off()

pdf(paste0(wd,"PhD_BPS53/plots/revision1/21_endothelium_scRNAseq/scRNAseq_DPT_meso_to_endo_celltypes_points_original.pdf"),
    width=5,height=2.7)
plot(x= tmp[ ,"DPT_meso_x1000"],
     y=as.numeric(as.factor(tmp$louvain_subclust)),
     col=alpha(all_colours[as.character(tmp$celltype)], 0.4), 
     ylab="",xlab="",pch=16,cex=1.2,ylim=c(0,5.2)
)
dev.off()




########==============================
# Plot ETS dynamics along pseudotime
########==============================


##Allantois trajectory

load(paste0(wd,"PhD_BPS53/data/revision1/21_endothelium_scRNAseq/snATACseq_endothelium_scRNAseq_ETS_allantois_dpt_varAbovepoint15_clusters_Normalised.rda"))
load(paste0(wd,"PhD_BPS53/data/revision1/21_endothelium_scRNAseq/snATACseq_endothelium_scRNAseq_ETS_allantois_varAbovepoint15_dpt_clusters_notNormalised.rda"))

clusters_al = clusters_notnorm

clust_al = data.frame(geneID = names(clusters_al$clusters),
           cluster = clusters_al$clusters,
           geneName = anSeq::getGeneName(names(clusters_al$clusters))[[2]])

clust_al[order(clust_al$cluster),]

library("tidyverse")

df0 = as.data.frame(clusters_al$predictions)

colnames(df0) = anSeq::getGeneName(colnames(df0))[[2]]
df0[,"sequence"] = seq(1,nrow(clusters_al$predictions))

df <- df0 %>%
  gather(key = "variable", value = "value", -sequence)
head(df)
clust_plot_ETS = list()

for (i in 1:length(unique(clust_al$cluster))){
  clust = unique(clust_al$cluster)[i]
  genesSel = as.character(clust_al$geneName[clust_al$cluster==clust])
  df_sub = df[as.character(df$variable) %in%as.character(genesSel),]
  plot = ggplot(df_sub, aes(x = sequence, y = value)) + 
    geom_line(aes(color = variable)) + ylim(-0.5,5) +# theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    ggtitle(paste("Allantois cluster",clust))  #
    #theme(axis.text = element_blank(), axis.title = element_blank(),
    #      axis.ticks = element_blank())
  clust_plot_ETS[[i]]=plot
  print(plot)
  
}
library(cowplot)
plot_fig = plot_grid(plotlist = clust_plot_ETS,ncol=3)
plot_fig

ggsave(plot_fig,
       file = paste0(wd,"PhD_BPS53/plots/revision1/21_endothelium_scRNAseq/scRNAseq_from_AL_to_endo_notnorm_clusters_as_in_gutPijuanSala_varAbovepoint15.pdf"),
       width = 8, height = 2)



##====
##Mixed mesoderm trajectory

load(paste0(wd,"PhD_BPS53/data/revision1/21_endothelium_scRNAseq/snATACseq_endothelium_scRNAseq_ETS_mesoderm_dpt_clusters_notNormalised_varAbovepoint15.rda"))
load(paste0(wd,"PhD_BPS53/data/revision1/21_endothelium_scRNAseq/snATACseq_endothelium_scRNAseq_ETS_mesoderm_dpt_clusters_Normalised_varAbovepoint15.rda"))

clusters_meso = clusters_notnorm
clust_meso = data.frame(geneID = names(clusters_meso$clusters),
                           cluster = clusters_meso$clusters,
                           geneName = anSeq::getGeneName(names(clusters_meso$clusters))[[2]])

clust_meso[order(clust_meso$cluster),]

library("tidyverse")

df0 = as.data.frame(clusters_meso$predictions)

colnames(df0) = anSeq::getGeneName(colnames(df0))[[2]]
df0[,"sequence"] = seq(1,nrow(df0))
                         
df <- df0 %>%
  gather(key = "variable", value = "value", -sequence)
head(df)


clust_plot_ETS = list()

for (i in 1:length(unique(clust_meso$cluster))){
  clust = unique(clust_meso$cluster)[i]
  genesSel = clust_meso$geneName[clust_meso$cluster==clust]
  df_sub = df[df$variable %in%as.character(genesSel),]
  plot = ggplot(df_sub, aes(x = sequence, y = value)) + 
    geom_line(aes(color = variable)) +
    ggtitle(paste("Mesoderm cluster",clust)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) + ylim(-0.5,5)
    #theme(axis.text = element_blank(), axis.title = element_blank(),
        #  axis.ticks = element_blank())
  clust_plot_ETS[[i]]=plot
  print(plot)
  
}

plot_fig = plot_grid(plotlist = clust_plot_ETS,ncol=3)
plot_fig
ggsave(plot_fig,
       file = paste0(wd,"PhD_BPS53/plots/revision1/21_endothelium_scRNAseq/scRNAseq_from_meso_to_endo_notnorm_clusters_as_in_gutPijuanSala_varAbovepoint15.pdf"),
       width = 8, height = 2)

