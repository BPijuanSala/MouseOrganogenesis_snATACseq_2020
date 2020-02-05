#!/usr/bin/perl


###############################################
##  Blanca Pijuan-Sala
##  Gottgens Lab
##  University of Cambridge, UK
##  12th December 2018
##  DESCRIPTION:
##  Program to reverse complement Barcode sequence as follows. The first
##  element of the fastq header corresponds to the unique barcode, consisting
##  fo the first 10 bp I1 + 12 bp I2. In HiSeq2500, we need to reverse reversecomplement
##  I2 to match the same cell in data produced from NextSeq500
##
##################################################


use strict;
use File::Basename;


if ( $#ARGV != 1 ) { print "\nPlease provide: <Input _density.txt> & <Genome(mm9/mm10/hg19/rn5)>\n"; exit; };

my $file = $ARGV[0]; chomp $file;
my $output = $ARGV[1]; chomp $output;

open FILE1, "$file";

open OUT, ">$output";

my $line_num = 0;
while (<FILE1>)
{
$line_num++;
my $line = $_;
chomp $line;

	#if ( $line =~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/ )
	if ($line =~ /^@.*\:.*/){

			my @spl = split(':', $line);
			my $origin_seq= reverse($spl[0]);

			my $origin_seq2 = substr($origin_seq, 0, -1);

			my $origin_seq3= reverse($origin_seq2);
			print "$spl[0]\n";
			print(length($origin_seq3));
			my $bc1 = substr($origin_seq3, 0, 10);
			print(length($bc1));

			my $bc2 = substr($origin_seq3, 10, 22);

			my $revcomp = reverse($bc2);

			$revcomp =~ tr/ATGCatgc/TACGtacg/;
			print "$revcomp\n\n";


			$spl[0] = "@" . $bc1 . $revcomp;

			#unshift @spl, '@';

			my $str = join ':', @spl;
			print (OUT "$str\n");
} else {
	print (OUT "$line\n");

}}
