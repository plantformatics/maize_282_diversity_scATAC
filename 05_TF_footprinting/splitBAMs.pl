#!/usr/bin/perl
use strict;
use warnings;

# load barcode information
print STDERR "...LOADING BARCODE METADATA\n";
my %gdata;
open F, $ARGV[0] or die;
while(<F>){
	chomp;
	if($_ =~ /^cellID/){
		next;
	}
	my @col = split("\t",$_);
	$col[9] =~ s/Goodman-Buckler//g;
	$gdata{$col[0]} = $col[9];
}
close F;

# load group info
print STDERR "...LOADING HAPLOTYPE PARTITIONING\n";
my @header;
my %acrdata;
open G, $ARGV[1] or die;
while(<G>){
	chomp;
	my @col = split("\t",$_);
	if($_ =~ /^282set/){
		@header = @col;
		next;
	}else{
		for (my $i = 1; $i < @col; $i++){
			my $idx = $i - 1;
			my @acr_coord = split("_",$col[0]);
			my $pos = join("-", @acr_coord[1..2]);
			$acrdata{$acr_coord[0]}{$pos}{$header[$idx]} = $col[$i];
		}
	}
}
close G;

# output files
my $celltype = $ARGV[2];
my $hap1out = $celltype . ".haplotype1.sam";
my $hap2out = $celltype . ".haplotype2.sam";
open (my $hap1, ">", $hap1out) or die;
open (my $hap2, ">", $hap2out) or die;

# iterate over stream
print STDERR "...STREAMING ALIGNMENTS\n";
my $its = 0;
my $pushed = 0;
open H, $ARGV[3] or die;
while(<H>){

	# verbose
	$its++;
	if(($its % 1000000) == 0){print STDERR "- iterated over $its records, pushing $pushed reads into haplogroups...\n";}

	# split
	chomp;
	my @col = split("\t",$_);

	# check if read matches a good BC
	if(! exists $gdata{$col[11]}){
		next;
	}

	# get read coordinates
	my $chr = $col[2];
	my $str = $col[3] - 1;
	my $end = $str + length($col[9]) + 1;
	
	# check if read falls in the correct interval
	my @sites = sort {$a cmp $b} keys %{$acrdata{$chr}};
	my $ol = 0;
	my $tracker = 0;
	for (my $i = 0; $i < @sites; $i++){
		my @pos = split("-",$sites[$i]);
		if($str <= $pos[1] && $end >= $pos[0]){
			$tracker++;
			$ol = $sites[$i];
			last;
			
		}
	}

	# if overlap an ACR
	if($tracker > 0 && $ol ne 0){
		my $geno = $gdata{$col[11]};
		my $happush = $acrdata{$chr}{$ol}{$geno};
		if($happush eq 1){
			$pushed++;
			print $hap1 "$_\n";
		}elsif($happush eq 2){
			$pushed++;
			print $hap2 "$_\n";
		}else{
			next;
		}
	}else{
		next;
	}
	
}
close H;

# close output
close $hap1;
close $hap2;
