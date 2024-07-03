#!/usr/bin/perl
use strict;
use warnings;

open F, $ARGV[0] or die;
while(<F>){
	chomp;
	if($_ =~ /^\@HD/){
		my @col = split("\t",$_);
		print "$col[0]\tVN:1.6\t$col[1]\n";
	}else{
		print "$_\n";
	}
}
close F;
