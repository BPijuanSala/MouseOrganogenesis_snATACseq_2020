#! /usr/bin/perl -w

use strict;

my $input = $ARGV[0];
chomp $input;

my $bp = $ARGV[1];
chomp $bp;

my $script_id = $ARGV[2];
chomp $script_id;

open ( IN, "$input" ) or die ( "Cannot open file $input \n" );

if ($script_id == 1)
{
open ( OUT_UP, ">50kb_upstream_of_peak.bed" );
open ( OUT_DOWN, ">50kb_downstream_of_peak.bed" );
}
elsif ($script_id == 2)
{
open ( OUT_UP, ">upstream_of_peak.bed" );
open ( OUT_DOWN, ">downstream_of_peak.bed" );
}

while(<IN>)
{
	my $line = $_;
	chomp $line;
	
	my @bits = split(/\t/, $line);
	my ($chr, $start, $end) = ($bits[0], $bits[1], $bits[2]);
	chomp($chr, $start, $end);

	my ($upstream_start, $upstream_end, $downstream_start, $downstream_end) = ("", "", "", "");

	$upstream_start = $start - $bp;
	if($upstream_start < 1){ $upstream_start = 1; }
	$upstream_end = $start;

	$downstream_start = $end;
	$downstream_end = $end + $bp;

	my $original_peak = $chr."_".$start."_".$end;

	print(OUT_UP "$chr\t$upstream_start\t$upstream_end\t$original_peak\n");
	print(OUT_DOWN "$chr\t$downstream_start\t$downstream_end\t$original_peak\n");
}	


close ( OUT_UP );
close ( OUT_DOWN );
close ( IN );
