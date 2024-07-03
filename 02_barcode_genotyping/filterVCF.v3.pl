#!/usr/bin/perl
use strict;
use warnings;
use POSIX;

# load genotypes
my %gen;
open G, $ARGV[1] or die;
while(<G>){
	chomp;
	$gen{$_} = 1;
}
close G;

# genotype indexes to keep
my @keep;

# iterate over VCF
open F, $ARGV[0];
while(<F>){
	chomp;
	if($_ =~ /^##/){
		print "$_\n";
		next;
	}elsif($_ =~ /#CHROM/){
		my @col = split("\t",$_);
		for (my $i = 9; $i < @col; $i++){
			if(exists $gen{$col[$i]}){
				push(@keep, $i);
			}
		}
		foreach(@col[0..7]){
			print "$_\t";
		}
		print "$col[8]";
		foreach(@keep){
			print "\t$col[$_]";
		}
		print "\n";
	}else{
		my @col = split("\t",$_);
		my $missing = 0;
		my $hetsites = 0;
		my $rhomsites = 0;
		my $ahomsites = 0;
		my $total = 0;
		my $alt = 0;
		for (my $j = 0; $j < @keep; $j++){
			my $idx = $keep[$j];
			if($j == -1){
				print STDERR " - place holder...\n";
			}else{
				if($col[$idx] =~ /0\|0/){
					$total = $total + 2;
					$rhomsites = $rhomsites + 2;
				}elsif($col[$idx] =~ /0\|1/){
					$total = $total + 2;
					$alt++;
					$hetsites = $hetsites + 2;
				}elsif($col[$idx] =~ /1\|1/){
					$total = $total + 2;
					$alt = $alt + 2;
					$ahomsites = $ahomsites + 2;
				}elsif($col[$idx] =~ /1\|0/){
					$total = $total + 2;
					$alt++;
					$hetsites = $hetsites + 2;
				}else{
					$missing++;
				}
			}
		}
		if($missing > 0){
			next;
		}elsif($total == 0){
			next;
		}elsif($hetsites > (int($total/2))){
			next;
		}elsif($rhomsites == 0 || $ahomsites == 0){
			next;
		}elsif($alt < 2){
			next;
#		}elsif($rhomsites <= 2 || $ahomsites <= 2){
		}else{
			my $af = sprintf("%.4f",$alt/$total);
			if($af > 0){
				$col[7] = "AF=$af;AC=$alt;AN=$total";
				foreach(@col[0..7]){
					print "$_\t";
				}
				print "$col[8]";
				foreach(@keep){
					print "\t$col[$_]";
				}
				print "\n";
			}
		}
	}
}
close F;
