#!/bin/bash
#SBATCH -A GOTTGENS-SL3-CPU
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p skylake
# set max wallclock time
#SBATCH --time=12:00:00
#SBATCH -e slurm.revComp.err

###############################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  April 2019
##  DESCRIPTION: Script to reverse complement I2 from HiSeq2500 fastq files.
###############################################


module load perl-5.20.0-gcc-5.4.0-4npvg5p

cd /path/to/directory/HiSeq2500_compl/fastq/

#Execute reverseComplement_barcodeheader_fastq.pl for each fastq file.
echo "Creating first line" >tmpFile_rev.txt
/path/to/directory/bin/reverseComplement_barcodeheader_fastq.pl /path/to/directory/HiSeq2500/fastq/HiSeq2500_all.demultiplexed.R1.fastq HiSeq2500_all_compl.demultiplexed.R1.fastq && echo "$i DONE" >> 'tmpFile_rev.txt' &
/path/to/directory/bin/reverseComplement_barcodeheader_fastq.pl /path/to/directory/HiSeq2500/fastq/HiSeq2500_all.demultiplexed.R2.fastq HiSeq2500_all_compl.demultiplexed.R2.fastq && echo "$i DONE" >> 'tmpFile_rev.txt' &
/path/to/directory/bin/reverseComplement_barcodeheader_fastq.pl /path/to/directory/HiSeq2500/fastq/HiSeq2500_large.demultiplexed.R1.fastq HiSeq2500_large_compl.demultiplexed.R1.fastq && echo "$i DONE" >> 'tmpFile_rev.txt' &
/path/to/directory/bin/reverseComplement_barcodeheader_fastq.pl /path/to/directory/HiSeq2500/fastq/HiSeq2500_large.demultiplexed.R2.fastq HiSeq2500_large_compl.demultiplexed.R2.fastq && echo "$i DONE" >> 'tmpFile_rev.txt' &
/path/to/directory/bin/reverseComplement_barcodeheader_fastq.pl /path/to/directory/HiSeq2500/fastq/HiSeq2500_small.demultiplexed.R1.fastq HiSeq2500_small_compl.demultiplexed.R1.fastq && echo "$i DONE" >> 'tmpFile_rev.txt' &
/path/to/directory/bin/reverseComplement_barcodeheader_fastq.pl /path/to/directory/HiSeq2500/fastq/HiSeq2500_small.demultiplexed.R2.fastq HiSeq2500_small_compl.demultiplexed.R2.fastq && echo "$i DONE" >> 'tmpFile_rev.txt' &


while (( "$( wc -l < tmpFile_rev.txt  | xargs)" < 7)); do
  wait
done

module unload perl-5.20.0-gcc-5.4.0-4npvg5p
