###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Generate Venn Diagram
###############################################################################


library(ggplot2)
wd <- "/path/to/directory/"


library(VennDiagram)

pdf(paste0(wd,"PhD_BPS53/plots/revision1/Heptad_peaks_VennDiagram.pdf"),
    width=5,height=5)
plot.new()

par(mar=c(4,4,4,4))
draw.pairwise.venn(2534,2370,750,fill=c("gray","#ff891c"),alpha=c(0.6,0.6))
dev.off()
