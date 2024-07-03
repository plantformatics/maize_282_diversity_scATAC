#!/usr/bin/perl
use strict;
use warnings;
#use Sort::Naturally;

my %flag = (
	0  => 'read_paired',
	1  => 'read_properly',
	2  => 'read_umapped',
	3  => 'mate_unmapped',
	4  => 'read_reverse',
	5  => 'mate_reverse',
	6  => 'first_pair',
	7  => 'second_pair',
	8  => 'secondary',
	9  => 'fail_qc',
	10 => 'duplicate',
	11 => 'sup_align'
);
	
my $counts = 0;
open F, "samtools view $ARGV[0] | " or die;
while(<F>){
	chomp;
	$counts++;
	if(($counts % 1000000) == 0){
		print STDERR " - iterated over $counts reads ... \n";
	}
	my @col = split("\t",$_);
	my $bin = sprintf ("%.12b", $col[1]);
	my @bins = split("",$bin);
	my @flags;

	# load flags
	for(my $i = 0; $i < @bins; $i++){
		if($bins[$i] == 1){
			push(@flags, $flag{$i});
		}
	}

	# get start and end pos of read
	my $chr = $col[2];
	my $pos1 = $col[3];
	my $cigar = parseCIGAR($col[5]);
	my %cig = %$cigar;
	my $netdif = 1;
	my @act = keys %cig;
	for (my $i = 0; $i < @act; $i++){
		my @action = split("_",$act[$i]);
		my $times = $cig{$act[$i]};
		if($action[1] eq 'M' || $action[1] eq 'S'){
			$netdif = $netdif + $times;
		}elsif($action[1] eq 'D'){
			$netdif = $netdif + $times;
		}
	}
	my $pos2 = $netdif + $pos1;
	if(grep {$_ =~ /read_reverse/} @flags){
		my $end = $pos2 - 4;
		my $start = $end - 1;
		print "$chr\t$start\t$end\t-\n";
	}else{
		my $start = $pos1 + 5;
		my $end = $start + 1;
		print "$chr\t$start\t$end\t+\n";
	}
}
close F;

# subroutines
sub extractcoord {

        # get read and position variables
        my ($hashref, $slot) = @_;
        my $readchr = $hashref->{'chr'};
        my $readpos = $hashref->{'pos'};
        my $readseq = $hashref->{$slot};
        my $readlen = length($readseq);

        # parse cigar
        my $readcigarref = parseCIGAR($hashref->{'cigar'});
        my %readcigar = %{$readcigarref};
        my @aln = keys %readcigar;

        # vars
        my %newseq;
        my @oldseq;
        if($slot eq 'seq'){
                @oldseq = split("",$readseq);
        }elsif($slot eq 'bq'){
                @oldseq = split("_",$readseq);
        }
        my $last = $readpos;
        my $tt = 0;
        my $penalty = 0;

        # iterate over CIGAR
        for (my $i = 0; $i < @aln; $i++){
                my @actiont = split("_",$aln[$i]);
                my $alntype = $actiont[1];
                if($alntype eq 'M'){
                        my $adj;
                        for (my $j = 0; $j < $readcigar{$aln[$i]}; $j++){
                                $adj = $j + $last;
                                $newseq{$adj} = $oldseq[$tt];
                                $tt++;
                        }
                        $last = $adj + 1;
                }elsif($alntype eq 'D'){
                        my $adj;
                        for (my $j = 0; $j < $readcigar{$aln[$i]}; $j++){
                                $adj = $j + $last;
                                $newseq{$adj} = '-';
                        }
                        $last = $adj + 1;
                }elsif($alntype eq 'I' || $alntype eq 'S'){
                        for (my $j = 0; $j < $readcigar{$aln[$i]}; $j++){
                                $tt++;
                        }
                }#elsif($alntype eq 'S'){
                #       my $adj;
                #       for (my $j = 0; $j < $readcigar{$aln[$i]}; $j++){
                                #$adj = $j + $last;
                #               $tt++;
                #       }
                        #$last = $adj + 1;
                #}
        }

        # return sequence
        my @newsites = sort {$a <=> $b} keys %newseq;
        my @outseq;
        for (my $k = 0; $k < @newsites; $k++){
                push(@outseq, $newseq{$newsites[$k]});
        }
        my $mergeseq = join("",@outseq);
        my $pos1 = $newsites[0];
        my $pos2 = $newsites[$#newsites];
        return($pos1, $pos2, \%newseq, $mergeseq);
}

sub parseCIGAR {
        my ($str) = @_;
        my @l = split(/(?<=\d)(?=\D)|(?<=\D)(?=\d)/, $str);
        my %hash;
        my $counter=0;
        for (my $i=0; $i < @l; $i+=2){
                my $action = $l[$i+1];
                my $times  = $l[$i];
                $counter++;
                my $id = join("_", $counter, $action);
                $hash{$id} = $times;
        }
        return(\%hash);
}

sub parseSAM {
        my ($line) = @_;
        my @col = split("\t",$line);
        my $strand;
        if($col[1] == 147){
                $strand = '-';
        }elsif($col[1] == 163){
                $strand = '+';
        }elsif($col[1] == 99){
                $strand = '+';
        }elsif($col[1] == 83){
                $strand = '-';
        }
        my @phred = split("", $col[10]);
        my @scores;
        foreach(@phred){
                my $phd = ord($_)-33;
                my $prob = (10**(-$phd/10));
                push(@scores, $prob);
        }
        my $score = join("_",@scores);
        my %hash = (
                'readname'      => $col[0],
                'flag'          => $col[1],
                'chr'           => $col[2],
                'pos'           => $col[3],
                'mq'            => $col[4],
                'cigar'         => $col[5],
                'pe_chr'        => $col[6],
                'pe_pos'        => $col[7],
                'insert'        => $col[8],
                'seq'           => $col[9],
                'bq'            => $score,
                'strand'        => $strand
        );
        foreach(@col[11..$#col]){
                if($_ =~ /CB/){
                        $hash{'cb'} = $_;
                }elsif($_ =~ /MD/){
                        $hash{'md'} = $_;
                }elsif(exists $hash{'cb'} && exists $hash{'md'}){
                        last;
                }
        }
        return(\%hash);
}

