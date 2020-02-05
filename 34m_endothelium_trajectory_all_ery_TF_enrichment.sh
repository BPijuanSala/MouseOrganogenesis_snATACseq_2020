#!/bin/bash
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -p cccp
# set max wallclock time
#SBATCH --time=2-00:00:00
#SBATCH -e slurm.HOMER_DPT.err

###############################################################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  October 2019
##  DESCRIPTION: Script to do TF motif enrichment analysis using HOMER
###############################################################################
cd /path/to/directory/sample_pooled_preprocess_revision1/18_endothelium_analysis/data/DPT_clusters_asGut

echo "files" >file.tmp

for i in {1..12}
do
  mkdir ./HOMER/'cluster_'$i && \
  /home/USSR/bp382/bin/Homer/bin/findMotifsGenome.pl 'snATACseq_endothelium_eryAl_DPT_cluster_'$i'.txt.bed' mm10 './HOMER/cluster_'$i'/' -size given && \
  echo "cluster $i done\n" >>file.tmp &

done


while (( "$( ls file.tmp | wc -l | xargs)" < 13)); do
  wait
done
