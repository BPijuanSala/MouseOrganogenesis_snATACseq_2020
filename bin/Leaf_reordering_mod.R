#####################################################################################
## file: Leaf_reordering.R
## author: Rebecca Hannah
## date: 15/11/2016
## description: It reorders the leaves of the hierarchical clustering so that they 
##              depending on their similarities
## Experiment: bin_scripts
#####################################################################################
leaf.reordering <- function(data, clust.method = "ward.D2",reorder="column",corMethod="spearman"){
  #reorder can be column or row.
  
  require(cba)
  if (reorder == "column"){
    
    mat1 <- as.matrix(data)
    
    cor.mat1<-cor(mat1,method=corMethod)
    dissim1<-(1-cor.mat1)/2
    d1<-as.dist(dissim1)
    
    #d <- dist(mat1)
    
    hc <- hclust(d1, method = clust.method, members=NULL)
    
    co <- order.optimal(d1, hc$merge)
    
    ho <- hc
    
    ho$merge <- co$merge
    
    ho$order <- co$order
    
    ho$labels
    return(ho)
    
  }
  if (reorder=="row"){
    
    mat2 <- t(mat1)
    
    
    cor.mat2<-cor(mat2,method=corMethod)
    dissim2<-(1-cor.mat2)/2
    d2<-as.dist(dissim2)
    
    #d2 <- dist(mat2)
    
    hc2 <- hclust(d2, method = clust.method, members=NULL)
    
    co2 <- order.optimal(d2, hc2$merge)
    
    ho2 <- hc2
    
    ho2$merge <- co2$merge
    
    ho2$order <- co2$order
    return(ho2)
  }
  if (reorder == "both"){
    
    mat1 <- as.matrix(data)
    
    cor.mat1<-cor(mat1,method=corMethod)
    dissim1<-(1-cor.mat1)/2
    d1<-as.dist(dissim1)
    
    #d <- dist(mat1)
    
    hc <- hclust(d1, method = clust.method, members=NULL)
    
    co <- order.optimal(d1, hc$merge)
    
    ho <- hc
    
    ho$merge <- co$merge
    
    ho$order <- co$order
    
    ho$labels
    
    
    mat2 <- t(mat1)
    
    
    cor.mat2<-cor(mat2,method=corMethod)
    dissim2<-(1-cor.mat2)/2
    d2<-as.dist(dissim2)
    
    #d2 <- dist(mat2)
    
    hc2 <- hclust(d2, method = clust.method, members=NULL)
    
    co2 <- order.optimal(d2, hc2$merge)
    
    ho2 <- hc2
    
    ho2$merge <- co2$merge
    
    ho2$order <- co2$order
    
    # require(roots)
    # require(gplots)
    # if(label.colors = NULL){
    
    #  plotgott <- heatmap.gottgens(mat1[genes,cells], Rowv = as.dendrogram(ho), Colv = as.dendrogram(ho2))
    #trace=c("none"), margins=c(14,14), key = TRUE, density.info=c("none"), sepwidth=c(0.01, 0.01), sepcolor="black", col = colorRampPalette(c("blue", "white", "red"))(100))
    
    #  } else {
    
    #   plotgott <- heatmap.gottgens(mat1[genes,cells], Rowv = as.dendrogram(ho), Colv = as.dendrogram(ho2),
    #                                ColSideColors = label.colors)
    # }
    return (list(cell.clust = ho,
                 gene.clust = ho2))
    
  }
}