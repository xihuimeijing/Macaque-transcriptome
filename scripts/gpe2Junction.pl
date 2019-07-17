#!/usr/bin/perl
##############################################################################################
# @version 1.0
# GPE to Junction file (chr	1base-lastExonEnd	1base-nextExonStart	name	strand)
#
# 							Sun Feb 19 11:05:22 CST 2012
#							Author: Jerry Chen			
##############################################################################################
use warnings;
use strict;
use Getopt::Long;
#use List::Util qw (sum min max);
use vars qw ($input $output $nb);
use vars qw ($fh1 $fh2);
GetOptions (
	"input|i:s"	=>	\$input,
	"output|o:s"	=>	\$output,
	"nobin|nb"	=>	\$nb,
	"help|h"        =>      sub{&usage;exit(0);}
);
unless ($input){$fh1 = \*STDIN;;}
else {open $fh1,$input;}
unless ($output) {$fh2 = \*STDOUT;}
else {open $fh2,">>$output";}



while(<$fh1>){
	chomp;
	next if (/^#/);
	my @field = split ("\t",$_);
	if (defined $nb){unshift @field,1;}
	my @starts = split (",",$field[9]);
	my @ends = split (",",$field[10]);
	
	if($field[3] eq '+'){
		for(my $i=0;$i<$#starts;$i++){
			my $lastEnd = $ends[$i];
			my $nextStart = $starts[$i+1]+1;
			my $state;
			if($field[6]<$lastEnd && $field[7]>$nextStart-1){
				$state = 'CDS-CDS';
			}elsif($field[6]>$nextStart-1){
				$state = '5UTR-5UTR';
			}elsif($field[7]<$lastEnd){
				$state = '3UTR-3UTR';
			}elsif($field[6]==$nextStart-1){
				$state = '5UTR-CDS';
			}elsif($field[7]==$lastEnd){
				$state = 'CDS-3UTR';
			}
			print $fh2 "$field[2]\t$lastEnd\t$nextStart\t";
			print $fh2 join('.',($field[1],$i,$state));
			print $fh2 "\t$field[3]\n";
		}		
	}else{
		for(my $i=0;$i<$#starts;$i++){
			my $lastEnd = $ends[$i];
			my $nextStart = $starts[$i+1]+1;
			my $state;
			if($field[6]<$lastEnd && $field[7]>$nextStart-1){
				$state = 'CDS-CDS';
			}elsif($field[6]>$nextStart-1){
				$state = '3UTR-3UTR';
			}elsif($field[7]<$lastEnd){
				$state = '5UTR-5UTR';
			}elsif($field[6]==$nextStart-1){
				$state = '3UTR-CDS';
			}elsif($field[7]==$lastEnd){
				$state = 'CDS-5UTR';
			}
			print $fh2 "$field[2]\t$lastEnd\t$nextStart\t";
			print $fh2 join('.',($field[1],$i,$state));
			print $fh2 "\t$field[3]\n";
		}
	}
}


sub usage{
print STDERR <<HELP
Usage: perl $0 -i [.GPE]  -o [output] -h
	--input|-i:s		input;GPE format;UCSC table
	--output|-o:s		output [Default:STDOUT]
	--nobin|-nb		gpe with no bin
	--help|-h		print this message
HELP
}
