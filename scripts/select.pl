#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use List::Util qw(max);
use vars qw ($input $reference $column1 $column2 $o);
use vars qw ($fh1 $fh2 $fh3);
GetOptions (
        "input|i:s"           =>      \$input,
	"reference|r:s"       =>      \$reference,
        "column1|c1:s"          =>      \$column1,
        "column2|c2:s"          =>      \$column2,
        "opposite|o:s"           =>     \$o,
	"help|h"              =>      sub{&usage;exit(0);}
);


open $fh1,"< $input";
open $fh2,"< $reference";
my ($lineA,$lineB,@array1,@array2,%hash);
$column1=$column1-1;
$column2=$column2-1;
while($lineA=<$fh1>)
    {
    chomp $lineA;
    @array1=split/\t/,$lineA;
    
    $hash{$array1[$column1]}=1;
    }
while($lineB=<$fh2>)
    {
    chomp $lineB;
    @array2=split/\t/,$lineB;
    if (defined $o)
        {          
        if (exists $hash{$array2[$column2]}){;}
        else {print "$lineB\n";}
        }
    else
    {
        if (exists $hash{$array2[$column2]}){print join ("\t",@array2);print "\n"}
    }
    }

sub usage{
print STDERR <<HELP
pick out the reference which contain the query
Usage: perl $0 -i [query]  -r [reference]  -c1 -c2 -o -h
        --input|-i:s             several column in which only column matters     
        --reference|-r:s         several column
        --column1|-c1:s          the column in query which matters
        --column2|-c2:s          the column that query lies in the refernce
        --opposite|-o            print the lines that does not exist in input
        --help|-h                print this message
HELP
}