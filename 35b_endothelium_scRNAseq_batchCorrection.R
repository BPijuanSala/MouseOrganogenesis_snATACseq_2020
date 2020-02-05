#############################################
#title: "10X embryo data - Characterising the haemato-endothelial landscape"
#author: "Blanca Pijuan-Sala"
#date: "09 September 2018"
#############################################

wd <- "/path/to/directory/scRNAseq_PijuanSala_2019/"
dir = "/path/to/directory/sample_pooled_preprocess_revision1/"
############################################
# FUNCTIONS
############################################

heatmapRedYelBlue <- c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090",
                       "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")

palette <- colorRampPalette(rev(heatmapRedYelBlue))
bluePal <- c("#BFBFBF","#6495ED","#000000")

####====================
# Load data
####====================

###Read metadata with the cell type annotation
meta = read.csv(paste0(wd,"data/PijuanSalaEtAl_SupplementaryTable3_revision.txt"), 
                header = T, stringsAsFactors = F, sep = "\t")


rownames(meta) <- meta$cell

###Load the counts of haemato-endothelial related clusters 

counts <- as.matrix(Matrix::readMM(file=paste0(dir,"21_endothelium_scRNAseq/data/countsData_endothelium_logged.mtx")))

colnames(counts) <- as.character(read.table(paste0(dir,"21_endothelium_scRNAseq/data/cells_endothelium.txt"),header=F)[,1])
rownames(counts) <-  read.table(paste0(dir,"21_endothelium_scRNAseq/data/genes_endothelium.txt"),header=F)[,1]

###Read metadata for the haemato-endothelial related clusters 

metaSub <- meta[colnames(counts),]
rownames(metaSub) <- metaSub$cell

####====================
# Batch correction
####====================

#functions
getHVGs = function(counts, min.mean = 1e-3, gene_df = genes){
  require(biomaRt)
  trend = scran::trendVar(counts,
                          loess.args = list(span = 0.05))
  decomp = scran::decomposeVar(counts, fit = trend)
  decomp = decomp[decomp$mean > min.mean,]
  #exclude sex genes
  xist = "ENSMUSG00000086503"
  mouse_ensembl = useMart("ensembl")
  mouse_ensembl = useDataset("mmusculus_gene_ensembl", mart = mouse_ensembl)
  # gene_map = getBM(attributes=c("ensembl_gene_id", "chromosome_name"), filters = "ensembl_gene_id", values = rownames(decomp), mart = mouse_ensembl)
  # ychr = gene_map[gene_map[,2] == "Y", 1]
  ychr = read.table(("/servers/bio-shares-4/gottgens/Blanca/Experiments/Rscripts/mmusculus_Y_genes.txt"), stringsAsFactors = FALSE)[,1]
  #other = c("tomato-td") #for the chimera
  decomp = decomp[!rownames(decomp) %in% c(xist, ychr),]
  decomp$FDR = p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < 0.05])
}

#ensure counts has columns names for the cells
#match timepoints,samples to the count table
#timepoint_order, sample_order should contain each sample/timepoint ONCE, in correct order for correction
doBatchCorrect = function(counts, timepoints, samples, timepoint_order, sample_order, npc = 50, BPPARAM = SerialParam()){
  require(scran)
  require(irlba)
  require(BiocParallel)
  #compute pca
  pca = prcomp_irlba(t(counts), n = npc)$x
  rownames(pca) = colnames(counts)
  if(length(unique(samples)) == 1){
    return(pca)
  }
  #create nested list
  pc_list = lapply(unique(timepoints), function(tp){
    sub_pc = pca[timepoints == tp, ]
    sub_samp = samples[timepoints == tp]
    list = lapply(unique(sub_samp), function(samp){
      sub_pc[sub_samp == samp, ]
    })
    names(list) = unique(sub_samp)
    return(list)
  })
  names(pc_list) = unique(timepoints)
  #arrange to match timepoint order
  pc_list = pc_list[order(match(names(pc_list), timepoint_order))]
  pc_list = lapply(pc_list, function(x){
    x[order(match(names(x), sample_order))]
  })
  #perform corrections within list elements (i.e. within stages)
  correct_list = lapply(pc_list, function(x){
    if(length(x) > 1){
      return(do.call(fastMNN, c(x, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected)
    } else {
      return(x[[1]])
    }
  })
  #perform correction over list
  if(length(correct_list)>1){
    correct = do.call(fastMNN, c(correct_list, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected
  } else {
    correct = correct_list[[1]]
  }
  correct = correct[match(colnames(counts), rownames(correct)),]
  return(correct)
}



hvgs = getHVGs(counts)


#get order: oldest to youngest; most cells to least cells
order_df = metaSub[!duplicated(metaSub$sample), c("stage", "sample")]
order_df$ncells = sapply(order_df$sample, function(x) sum(metaSub$sample == x))
order_df$stage = factor(order_df$stage,
                        levels = rev(c("E8.5",
                                       "E8.25",
                                       "E8.0",
                                       "E7.75",
                                       "E7.5",
                                       "E7.25",
                                       "mixed_gastrulation",
                                       "E7.0",
                                       "E6.75",
                                       "E6.5")))
order_df = order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
order_df$stage = as.character(order_df$stage)



cat("Batch correction...")

all_correct = doBatchCorrect(counts = counts[rownames(counts) %in% hvgs,],
                             timepoints = metaSub$stage,
                             samples = metaSub$sample,
                             timepoint_order = order_df$stage,
                             sample_order = order_df$sample,
                             npc = 50)
dim(all_correct)
Matrix::writeMM(Matrix(all_correct, sparse = TRUE),file=paste0(dir,"21_endothelium_scRNAseq/data/counts_endothelium_PCAcorrected.mtx"))


