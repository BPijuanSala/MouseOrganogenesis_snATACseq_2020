#! /usr/bin/perl -w

use strict;
use autodie;

print "\n\nIs module bedtools/2.20.1 loaded? (y/n):";
my $answer = <STDIN>;
chomp($answer);
unless ($answer eq "y") {print "\n\nPlease load using \"module load bedtools/2.20.1\"\n\n"; exit;}

## gets genes whose PROMOTERS directly overlap peaks
## then genes DIRECTLY overlapping peaks
## then the closest gene either side WITHIN 50KB (if there is one!)

my $peak_in = $ARGV[0];

my $output = $ARGV[1];

my $input = $ARGV[2];

my $source = $ARGV[3];

my $prom_db = $ARGV[4];


if (@ARGV !=  5){ print("\nPlease provide as arguments; input file, output file, genome selection (mm10/hg19/), source selection (Ensembl = E, or UCSC = U - currently not available), promoter selection ('1000bp upstream of TSS' = T, or MPromDB = M) .\n\n"); exit; }

chomp ($peak_in, $output, $input, $source, $prom_db);

my $command = "";

if ($source eq "E")
{
	print "\nRunning using Ensembl refGene coordinates...\n";
	$command = "/path/to/file/bin/peaktogenes_dependencies/peaktogenes_Ensembl_v2.pl $peak_in $output $input $prom_db";
}
else
{
	print "\nSource selection not recognised or not available. Please specify "E".\n\n";
	exit;
}

system($command);

exit;
