#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
my ($infile,$a,$n,$max,$min,$outformat,$outnumber,%output,$e);

GetOptions(
	   'f|infile=s'       => \$infile,
	   'value=i'		=>\$a,
	   'n=s'		=>\$n,
	   'max=i'		=>\$max,
	   'min=i'		=>\$min,
	   'equal=s'		=>\$e,
	   'outbin=s'		=>\$outformat,
	   'outnumber=s'	=>\$outnumber,
           'h|help'          => sub{usage()}
	  ) || usage();


my $allbin=$max-$min; my $bin=$allbin/$n; my $all=0;
for (my $i=1;$i<=$n;$i++){$output{$i}=0;}
open INDEX,"$infile"or die "Can't open file $infile:$!";
while (<INDEX>){
    chomp;
    my @line=split (/\t/,$_); my $value=$line[$a-1]; $all++;
    for (my $i=1;$i<=$n;$i++){
	my $a=$min+($i-1)*$bin; my $b=$a+$bin;
	if ($e eq "max"){ if ($value>$a && $value<=$b){ $output{$i}+=1 } }
	if ($e eq "min"){ if ($value>=$a && $value<$b){ $output{$i}+=1 }  }
	 #if ($value>$a && $value<=$b){ $output{$i}+=1 }
    }
}    
if ($outnumber eq "percent"){for (my $i=1;$i<=$n;$i++){ $output{$i}=$output{$i}/$all } }
for (my $i=1;$i<=$n;$i++){
    my $a=$min+($i-1)*$bin; my $b=$a+$bin;
    if ($outformat eq "max") { say $b,"\t",$output{$i} }
    if ($outformat eq "min") { say $a,"\t",$output{$i} }
    if ($outformat eq "min~max") { say $a,"~",$b,"\t",$output{$i} }
    if ($outformat eq "number") { say $i,"\t",$output{$i} }
}


sub usage{
print <<HELP;
Usage:  perl $0 -f file -n 10 -max 1 -min 0 >file
Author:	Shuxian Li 2017-09-28
Options:
	-f	File	the file contain value
	-value	2 	the lie number of the value
	-n	10      the number of bin
	-max	1	the max of bin
	-min	0	the min of bin
	-equal	max or min 			  the boundary equal max or min 
	-outbin max or min or min~max or number	  the output bin format
	-outnumber number or percent	     the output value if number or percent
	-h    --help	Print this help information.
Output: bin1 percent
	bin2 percent ...
HELP
    exit(-1);
}
