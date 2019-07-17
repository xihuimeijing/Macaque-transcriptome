#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
my ($infile,$bin,$value,@line,$n,$max,$min,$outformat,$outnumber,%output);

GetOptions(
	   'f|infile=s'       => \$infile,
           'h|help'          => sub{usage()}
	  ) || usage();

open INDEX,"$infile"or die "Can't open file $infile:$!";
while (<INDEX>){
    chomp;
    @line=split (/\t/,$_);
    $bin=$line[0];$value=$line[1];
    if (exists $output{$bin}){
	$output{$bin}=$output{$bin}.$value.",";
    }else{$output{$bin}=$value.",";}
}

while (($bin,$value)=each %output){
    say $bin,"\t",$value;
}




sub usage{
print <<HELP;
Usage:  perl $0 -f file  >file
Author:	Shuxian Li 2018-03-24
Options:
	-f	File	
	-h    --help	Print this help information.
Input File: bin1 value1
	    bin1 value2
	    bin2 value3
	    bin3 value4
	    bin3 value5
Output: bin1 value1,value2
	bin2 value3
	bin3 value4,value5
HELP
    exit(-1);
}
