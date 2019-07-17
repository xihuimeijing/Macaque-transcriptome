#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
my($inFile,$column);
$column=1;
my $opt=GetOptions(
                        'i|in=s'         => \$inFile,
                        'c|col=i'        => \$column,
                        'h|help'        => sub{&usage;exit(-1);}
                  );
my ($IN,%number,%hash);
if(defined $inFile){
	open $IN,"$inFile" or die "Can't open file $inFile:#!";
}else{
	$IN=\*STDIN;
}
while(<$IN>){
	chomp;
	my @split=split;
	if(! defined $number{$split[$column-1]}){
		$number{$split[$column-1]}=1;
		$hash{$split[$column-1]}=$_;
	}else{
		$number{$split[$column-1]}+=1;
		$hash{$split[$column-1]}=$_;
	}
}
foreach my $key(sort keys %number){
	if($number{$key}==1){
		say $hash{$key};
	}
}
sub usage{
print STDERR <<HELP 
Usage:  perl $0 -i <IN.tsv> -c [1] >OUT.tsv
        Output unique result based on the values of specific column. That is remove rows with redundant values.
        'i|in'      FILE   The tab-delimit input file(If not provided, it can be read from STDIN).
        'c|col'     INT    The column number with redundant values to remove.(default:1)
        'help|h'        print this help message    
HELP
}