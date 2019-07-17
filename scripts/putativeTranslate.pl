#!/usr/bin/perl
################################################################################
#	@version 1.0
#		function: translate a sequence in 6-frames and choose the longsst ORF
#		input: .fasta
#
#						Author: Jerry Chen
#						Sat Nov 24 20:25:52 CST 2012
###############################################################################
use strict;
use warnings;
use lib "/rd1/user/chenjy/perl_lib/lib/perl5/";
use Bio::Tools::CodonTable;
use Getopt::Long qw(:config no_ignore_case bundling);
use List::Util qw (max);
use vars qw ($fh $fasta $frame $atg);

($frame,$atg) = (3);
my ($outnuc,$outpro,$outinfo) = ("long_nuc.rst","long_pro.rst","putativeCDS");
open OUT1,">>$outnuc";
open OUT2,">>$outpro";
open OUT3, ">>$outinfo";

GetOptions(
	"fasta|i:s"		=>	\$fasta,
	"frame|f:i"		=>	\$frame,
	"atg|a"			=>	\$atg,
	"help|h"		=>	sub{&usage;exit(0);}
);

# intialize the codon table
my $codonTable = Bio::Tools::CodonTable->new();$codonTable->id(1);	# standard codon table

unless($fasta){$fh =\*STDIN;}
else{open $fh,$fasta;}
if ($frame != 3 &&  $frame != 6){print STDERR "only 3 or 6 are accepted by --frame|-f\n";exit(0);}

my ($id,$nextid,$seq)=("","","");
my (@aa,@nuc);


while(<$fh>){
	chomp;
	if (/^>(.*)$/){$id = $1;}
	else{$id = $nextid;$seq = $_;}
	($seq,$nextid) = &retrive_block(my $newline = <$fh>);
# initialize the variables! 
	@aa = ();
	@nuc =();
# six-frame translation
	if($frame >= 3){
		my ($tmp_aa,$tmp_nuc);
		($tmp_aa,$tmp_nuc) = &longestORF($seq);
		push @aa,[$tmp_aa,1];
		push @nuc,[$tmp_nuc,1];
		($tmp_aa,$tmp_nuc) = &longestORF(my $seq1 = substr $seq,1,length($seq)); 
		push @aa,[$tmp_aa,2];
		push @nuc,[$tmp_nuc,2];
		($tmp_aa,$tmp_nuc) = &longestORF(my $seq2 = substr $seq,2,length($seq));  
		push @aa,[$tmp_aa,3];
		push @nuc,[$tmp_nuc,3];
		if ($frame == 6){
			my $revseq = join('',reverse(split('',$seq)));
			$revseq =~ tr /AGCTNagctn/TCGANTCGAN/;
			($tmp_aa,$tmp_nuc) = &longestORF($revseq);
			push @aa,[$tmp_aa,-1];
			push @nuc,[$tmp_nuc,-1];
			($tmp_aa,$tmp_nuc) = &longestORF(my $seq1 = substr $revseq,1,length($revseq)); 
			push @aa,[$tmp_aa,-2];
			push @nuc,[$tmp_nuc,-2];
			($tmp_aa,$tmp_nuc) = &longestORF(my $seq2 = substr $revseq,2,length($revseq));  
			push @aa,[$tmp_aa,-3];
			push @nuc,[$tmp_nuc,-3];
		}
	}
	my $index = max_index(my @length = map{length($_->[0])}@aa);
	
	my ($start1,$end1,$strand);
	$strand = '+';
	if ($aa[$index][1] < 0){
		($seq = join('',reverse(split('',$seq)))) =~ tr /AGCTNagctn/TCGANTCGAN/;
		$strand = '-';
	}
	if($seq =~ /($nuc[$index][0])/){
		$start1 = $-[0] + 1;
		$end1 = $+[0];
	}
	
	unless(defined $start1){$start1 = '-';$end1 = '-';}
	
	print OUT1 ">$id,$start1,$end1,$strand\n";
	print OUT1 "$nuc[$index][0]\n";
	print OUT2 ">$id\n$aa[$index][0]\n";
	print OUT3 "$id\t$start1\t$end1\t$strand\n";
}


sub retrive_block{
	if(my $newline = $_[0]){
		chomp ($newline);
		if ($newline =~ /^>(.*)$/){
			$nextid = $1;
		}else{
			$seq .= $newline;
			&retrive_block(my $newline = <$fh>);
		}
	}
	return(($seq,$nextid));
}

# input:	nucleotide sequence
# output:	1)the AA sequence
# 		2)the nucleotide sequence
# 		of the longest ORF (if there are more than one longest ORFs, only report the first one in this version)
sub longestORF{
	my $nucleotide_seq = $_[0];
	my $aa_seq = $codonTable->translate($nucleotide_seq);
	my @aaSequence; my @starts; my @ends;

	if ($atg){
		@aaSequence = $aa_seq =~ /(M[^*]*\*)/g;
		for(;$aa_seq =~ /(M[^*]*\*)/g;){
			push @starts,$-[0];
			push @ends,$+[0];
		}
		
	}else{
		@aaSequence = $aa_seq =~ /([^*]+\*)/g;
		for(;$aa_seq =~ /([^*]+\*)/g;){
			push @starts,$-[0];
			push @ends,$+[0];
		}
	}
	if (exists $aaSequence[0]){
		my $index = &max_index(my @length = map {$ends[$_]-$starts[$_]}(0..$#starts));
		return($aaSequence[$index],(substr $nucleotide_seq,3*$starts[$index],3*$length[$index]));
	}else{
		return(('-','---'));
	}
}

# OK! the first index for the max number!!!!!!!!!!!!!!!!!!!
sub max_index{
	my $max = max(@_);
	my @index = map{$_}grep{$_[$_]==$max}(0..$#_);
	return($index[0]);
}

sub usage {
print STDERR <<HELP
perl $0 --fasta|-i [input.fasta] --frame|-f [3/6] --atg|-a --help|-h
Options:
	--fasta|-i	input file with fasta format 
	--frame|-f	3-frame or 6-frame translation [default:3]
	--atg|-a	only aa sequence with M is considered
	--help|-h	print this help information
HELP
}
