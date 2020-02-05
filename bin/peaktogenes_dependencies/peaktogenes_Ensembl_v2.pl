#! /usr/bin/perl -w

use strict;

## modified peaktogenes1 considering strand:
## first gets genes with PROMOTERS overlapping peaks
## then genes DIRECTLY overlapping peaks
## then the 2 closest genes with WITHIN 50KB (either side) - if there is one!

my $peak_in = $ARGV[0];

my $output = $ARGV[1];

my $input = $ARGV[2];

my $prom_db = $ARGV[3];


if (@ARGV !=  4){ print("\nPlease provide as arguments; input file, output file, genome selection (mm9/mm10/hg18/hg19/rn5) - only mm10 available, promoter selection ('1000bp upstream of TSS' = T, or MPromDB = M) - only T available .\n\n"); exit; }

chomp ($peak_in, $output, $input, $prom_db);

my $script_id = 1;


##### Get genes info #####

my $gene_in = "";
my %gene = ();


### Get Ensembl genes ###

print "\nGetting ensembl genes...\n";

if ($input eq "mm10")
{
$gene_in = "/path/to/file/bin/mm10_geneEnsembl_sorted.txt";
}
else
{
print "\nGenome selection not recognised or not available\n";
exit;
}


## Prepare genes ##

print "\nPreparing genes...\n";

chomp $gene_in;

open(IN1, $gene_in) || die("Cannot open \$gene_in: $!.\n");

while(<IN1>)
{
my $templine2 = $_;
chomp $templine2;

	if ( $templine2 =~ /^#/ ) { next; }

	my @content = split(/\t/, $templine2);
	my $chr = $content[0];
	my $start = $content[1];
	my $end = $content[2];
	my $ensembl_ID = $content[4];
	my $name = $content[5];
	my $strand = $content[3];
	my $key = $chr."\t".$start."\t".$end;
	my $new_key;

	if ($gene{$key})
	{
	$new_key = $key."\t2";
	@{$gene{$new_key}} =  ($ensembl_ID, $name, $strand);
	}
	else
	{
	@{$gene{$key}} =  ($ensembl_ID, $name, $strand);
	}

}#close while

close(IN1);

my $prom_db_file = "";


##### Get promoter coordinates #####

print "\nGetting promoter coordinates...\n";

if ($input eq "mm10")
{
	if ($prom_db eq "T")
	{
		##THIS IS THE ONE THAT I MODIFIED!!!!
	$prom_db_file = "/path/to/file/bin/mm10_geneEnsembl_sorted_gene1kbup.txt";
	}
}



##### Go through each peak and compare to genes #####

#print "\nmarker1\n";


##### First check whether the peaks overlap promoter regions #####

#print "\nmarker2\n";

print "\nChecking peak overlaps with promoters...\n";

my $command = "intersectBed -a $peak_in -b $prom_db_file -wa -wb > promoter_overlaps_temp.bed";
system($command);

$command = "intersectBed -a $peak_in -b $prom_db_file -wa -wb > promoter_overlaps_final.bed";
system($command);


$command = "intersectBed -a $peak_in -b $prom_db_file -v > unmapped_peaks_temp.bed";
system($command);

my $unmapped = "unmapped_peaks_temp.bed";


## If using MPromDB, modify promoter overlaps - match gene names to ensembl IDs

### For non-overlapping peaks, get those directly overlapping genes and peaks not directly overlapping a gene

print "\nChecking direct peak overlaps with genes...\n";

$command = "intersectBed -a $unmapped -b $gene_in -wa -wb > gene_overlaps_temp.txt";
system($command);

$command = "intersectBed -a $unmapped -b $gene_in -wa -wb > gene_overlaps_final.txt";
system($command);

$command = "intersectBed -a $unmapped -b $gene_in -v > unmapped_peaks2_temp.bed";
system($command);

$unmapped = "unmapped_peaks2_temp.bed";



### Process remaining unmapped peaks ###

#print "\nmarker3\n";

print "\nChecking for genes within 50kb around peaks...\n";

#get 50kb up and downstream of each peak
$command = "perl /path/to/file/bin/peaktogenes_dependencies/get_extended_coordinates_dep.pl $unmapped 50000 $script_id";
system($command);

#check for gene overlaps within 100kb
$command = "intersectBed -a 50kb_upstream_of_peak.bed -b $gene_in -wa -wb > 50kb_up_overlaps_temp.bed";
system($command);

$command = "intersectBed -a 50kb_upstream_of_peak.bed -b $gene_in -v > unmapped_peaks3_temp.bed";
system($command);

$command = "intersectBed -a 50kb_downstream_of_peak.bed -b $gene_in -wa -wb > 50kb_down_overlaps_temp.bed";
system($command);

$command = "intersectBed -a 50kb_downstream_of_peak.bed -b $gene_in -v > unmapped_peaks4_temp.bed";
system($command);

#File remaining unmapped peaks

#print "\nmarker4\n";

print "\nWriting unmapped peaks file...\n";

$command = "cut -f4 unmapped_peaks3_temp.bed > unmapped_peaks5_temp.bed";
system($command);

$command = "cut -f4 unmapped_peaks4_temp.bed > unmapped_peaks6_temp.bed";
system($command);

$command = "cat unmapped_peaks5_temp.bed unmapped_peaks6_temp.bed | sort | uniq -d | sed 's/_/\t/g' > UNMAPPED_PEAKS.txt";
system($command);



### Find NEAREST up and downstream gene for each peak ###

#Get all peaks
#$command = "cat 50kb_up_overlaps_temp.bed 50kb_down_overlaps_temp.bed > temp.txt";
#system($command);
#$command = "cut -f4 temp.txt | sort | uniq > temp2.txt";
#system($command);

#open (IN_TEMP, "temp2.txt") or die;
#my @peaks = <IN_TEMP>;
#close(IN_TEMP);



### Find nearest gene WITHIN 50KB EITHER SIDE ###

print "\nFinding nearest gene within 50kb either side...\n";

#Process up and downstream gene overlaps
open (IN_UP, "50kb_up_overlaps_temp.bed") or die;
open (IN_DOWN, "50kb_down_overlaps_temp.bed") or die;

my %up_hash = ();
my %down_hash = ();
my $key = "";

#For each upstream overlap...
while(<IN_UP>)
{
my $distance = "";
my $line1 = $_;
chomp $line1;
my ($peak_chr_50kb, $peak_start_50kb, $peak_end_50kb, $original_peak, $gene_chr, $gene_start, $gene_end, $strand, $ensembl, $gene) = split(/\t/, $line1);

$key = $original_peak;
my ($peak_chr, $peak_start, $peak_end) = split(/\_/, $original_peak);

#print "\noriginalkey:$key\n";

#Get distance between peak start and gene start/end (depending on strand)
if ($strand eq "1")
{
$distance = $peak_start - $gene_start;
}
else
{
$distance = $peak_start - $gene_end;
}

if ($up_hash{$key})
{
my $value1 = $up_hash{$key};
my @values = split(/\t/, $value1);
my $value = shift(@values);
#print "\nlastvalue = $value\n";
	if ($distance < $value)
	{
	$up_hash{$key} = $distance."\t".$ensembl."\t".$gene."\t".$strand;
	}
#print "\nreplacedvalue = $up_hash{$key}\n";
}
else
{
$up_hash{$key} = $distance."\t".$ensembl."\t".$gene."\t".$strand;
#print "\nnewvalue = $up_hash{$key}\n";
}

}

#For each downstream overlap
while(<IN_DOWN>)
{
my $distance = "";
my $line2 = $_;
chomp $line2;
my ($peak_chr_50kb, $peak_start_50kb, $peak_end_50kb, $original_peak, $gene_chr, $gene_start, $gene_end, $strand, $ensembl, $gene) = split(/\t/, $line2);

$key = $original_peak;
my ($peak_chr, $peak_start, $peak_end) = split(/\_/, $original_peak);


#Get distance between peak start and gene start/end (depending on strand)
if ($strand eq "1")
{
$distance = $gene_start - $peak_end;
}
else
{
$distance = $gene_end - $peak_end;
}

if ($down_hash{$key})
{
my $value1 = $down_hash{$key};
my @values = split(/\t/, $value1);
my $value = shift(@values);
	if ($distance < $value)
	{
	$down_hash{$key} = $distance."\t".$ensembl."\t".$gene."\t".$strand;
	}
}
else
{
$down_hash{$key} = $distance."\t".$ensembl."\t".$gene."\t".$strand;
}

}

#Each hash should now contain only the nearest overlap for each peak - print


open (OUT_UP, ">50kb_up_overlaps2_temp.bed") or die;
open (OUT_DOWN, ">50kb_down_overlaps2_temp.bed") or die;

foreach my $key1 (keys %up_hash)
{
my ($dist, $ensembl_id, $gene_name, $strand) = split(/\t/, $up_hash{$key1});
#if ($dist > 50000){print(UNMAPPED_UP "$key1\n"); next;} #dist shouldnt be >50kb as we only searched within 50kb!!!
print (OUT_UP "$key1\t$ensembl_id\t$gene_name\t$strand\t$dist\n");

#print "\nprintkey:$key1\n";
#print "\nprintvalue:$up_hash{$key}\n";
}

foreach my $key2 (keys %down_hash)
{
my ($dist, $ensembl_id, $gene_name, $strand) = split(/\t/, $down_hash{$key2});
#if ($dist > 50000){print(UNMAPPED_DOWN "$key2\n"); next;} #dist shouldnt be >50kb as we only searched within 50kb!!!
print (OUT_DOWN "$key2\t$ensembl_id\t$gene_name\t$strand\t$dist\n");
}


#close(UNMAPPED);

#$command = "sort UNMAPPED_PEAKS_temp.bed | uniq > UNMAPPED_PEAKS.bed";
#system($command);

#Replace underscores in 50kb overlaps with tabs

$command = "sed 's/_/\t/g' 50kb_up_overlaps2_temp.bed > 50kb_up_overlaps2_temp2.bed";
system($command);

$command = "sed 's/_/\t/g' 50kb_up_overlaps2_temp.bed > 50kb_up_overlaps2_final.bed";
system($command);
$command = "sed 's/_/\t/g' 50kb_down_overlaps2_temp.bed > 50kb_down_overlaps2_temp2.bed";
system($command);

#ADDED BPS
$command = "sed 's/_/\t/g' 50kb_down_overlaps2_temp.bed > 50kb_down_overlaps_final.bed";
system($command);


print "\nPreparing output files...\n";

#Concatenate all matches into one file
$command = "cut -f8-9 promoter_overlaps_temp.bed > genes_temp.txt";
system($command);

$command = "cut -f1-3 promoter_overlaps_temp.bed > coords_temp.txt";
system($command);

$command = "cut -f7 promoter_overlaps_temp.bed > strand_temp.txt";
system($command);


$command='lines=$(cat strand_temp.txt | wc -l) &&  seq $lines | sed "c TSS" >distance_temp.txt';
system($command);
$command = "paste coords_temp.txt genes_temp.txt strand_temp.txt distance_temp.txt > promoter_overlaps2_temp.bed";
system($command);





$command = "cut -f8-9 gene_overlaps_temp.txt > genes_temp.txt";
system($command);

$command = "cut -f1-3 gene_overlaps_temp.txt > coords_temp.txt";
system($command);

$command = "cut -f7 gene_overlaps_temp.txt > strand_temp.txt";
system($command);



$command='lines=$(cat strand_temp.txt | wc -l) &&  seq $lines | sed "c gene" >distance_temp.txt';
system($command);

$command = "paste coords_temp.txt genes_temp.txt strand_temp.txt distance_temp.txt > gene_overlaps2_temp.bed";
system($command);

$command = "cat promoter_overlaps2_temp.bed gene_overlaps2_temp.bed 50kb_up_overlaps2_temp2.bed 50kb_down_overlaps2_temp2.bed > $output";
system($command);


open(IN_LAST, "$output") or die;
open(OUT_LAST, ">temp3.txt") or die;

while(<IN_LAST>)
{
my $line = $_;
chomp ($line);
$line =~ s/\_/\t/g;
chomp ($line);
print (OUT_LAST "$line\n");
}

$command = "mv temp3.txt $output";
system($command);


## Extract unique genes to produce second output file:

print "\nExtracting unique gene lists...\n";

## Gene name list
my $cmd = "cut -f5 $output > temp4.txt";
system($cmd);
my @filename = split(/\./, $output);
pop(@filename);
my $new_filename = join(".", @filename);
chomp $new_filename;
my $output2 = $new_filename."_names_unique.txt";
$cmd = "tr '[:lower:]' '[:upper:]' < temp4.txt | sort | uniq -i > $output2";
system($cmd);

## Ensembl ID list
$cmd = "cut -f4 $output > temp4.txt";
system($cmd);
@filename = split(/\./, $output);
pop(@filename);
$new_filename = join(".", @filename);
chomp $new_filename;
$output2 = $new_filename."_ensemblIDs_unique.txt";
$cmd = "tr '[:lower:]' '[:upper:]' < temp4.txt | sort | uniq -i > $output2";
system($cmd);


## Remove duplicates entries from final output
$command = "sort $output | uniq > temp_final.txt";
system($command);
$command = "mv temp_final.txt $output";
system($command);


## Remove temp files

#$cmd = "rm *temp*";
#system($cmd);

$cmd = "rm 50kb_downstream_of_peak.bed";
system($cmd);

$cmd = "rm 50kb_upstream_of_peak.bed";
system($cmd);

print "\nFINISHED.\n\n";

exit;
