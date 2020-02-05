###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Endothelium - plot enriched motifs
###############################################################################


library(ggplot2)
wd <- "/path/to/file/"


df0 = xlsx::read.xlsx(paste0(wd,"PhD_BPS53/data/revision1/18_endothelium_analysis/cluster_DPT/patterns_ATACseq.xlsx"),
                sheetIndex=1)
df = df0
df[] <- lapply(df, gsub, pattern='\t', replacement='')

colnames(df)[1] = "TF motif"
rownames(df) = df$`TF motif`
df = df[,-1]
colnames(df) = paste("pattern",seq(1:12))

df = as.matrix(df)
df = apply(df, 2, as.numeric)
rownames(df) = df0$NA.

heatmapRedYelBlue <- c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090",
                       "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")

palette <- colorRampPalette(rev(heatmapRedYelBlue))
pdf(paste0(wd,"PhD_BPS53/plots/revision1/18_endothelium_analysis/TF_motif_enrichment_patterns_heatmap.pdf"),
    width=5,height=7)
gplots::heatmap.2(-log(df), trace="none", 
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
