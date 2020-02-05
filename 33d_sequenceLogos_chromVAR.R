###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to plot sequence logos
###############################################################################

wd <- "/path/to/directory/"



library(chromVARmotifs)
library(PWMEnrich)


data("mouse_pwms_v2")

mouse_pfms_v3 = mouse_pwms_v2
for (i in names(mouse_pwms_v2@listData)){
  mouse_pfms_v3@listData[[i]] = 0.25*(exp(as.matrix(mouse_pwms_v2@listData[[i]])))
  }

mouse_pfms_v4 =mouse_pfms_v3
for (i in names(mouse_pfms_v3@listData)){
  mouse_pfms_v4@listData[[i]] =round(mouse_pfms_v3@listData[[i]]*100,0)
  mode(mouse_pfms_v4@listData[[i]]) <- "integer"
}

library(seqLogo)

for (i in names(mouse_pfms_v3@listData)){
  
  png(paste0(wd,"PhD_BPS53/plots/TF_motif_logos/",i,".png"),
      width = 7, height = 3, units = 'in', res = 300)
  seqLogo(as.matrix(mouse_pfms_v3@listData[[i]]),ic.scale = T)
  dev.off()
}

for (i in names(mouse_pfms_v3@listData)){
  
  pdf(paste0(wd,"PhD_BPS53/plots/TF_motif_logos/",i,".pdf"),
      width = 7, height = 4, units = 'in', res = 300)
  seqLogo(as.matrix(mouse_pfms_v3@listData[[i]]),ic.scale = T)
  dev.off()
}