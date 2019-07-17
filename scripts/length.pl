#!perl
#use strict;

use List::Util qw (sum min max);

use Getopt::Long;
use vars qw ($input $bin);
use vars qw ($fh1);
GetOptions (
         "input|i:s"                  =>      \$input,
         "bin|b"                      =>      \$bin,
         "help|h"                     =>      sub{&usage;exit(0);}
);

if (! open $fh1,"< $input"){die "cannot open $input!";}

$n=0;$sum=0;
while($line=<$fh1>)
	{
	$length=0;
	chomp $line;
	@array=split/\t/,$line;
	shift @array if defined $bin;
	@start=split/,/,$array[8];
	@end=split/,/,$array[9];
	for ($i=0;$i<$array[7];$i++)
		{
		if (($start[$i]<=$array[5])&&($end[$i]>=$array[6])){$startexon=$i;$endexon=$i;$length=$array[6]-$array[5];last;}
		if (($start[$i]<=$array[5])&&($end[$i]>=$array[5])){$startexon=$i;$length=$length+$end[$i]-$array[5];}
		if (($start[$i]<=$array[6])&&($end[$i]>=$array[6])){$endexon=$i;$length=$length+$array[6]-$start[$i];last;}
		}
	for ($i=$startexon+1;$i<$endexon;$i++)
		{
		$length=$length+$end[$i]-$start[$i];
		}
        for ($i=0;$i<$array[7];$i++)
                {
                $genelength=sum(@end)-sum(@start);
                }
	print "$array[0]\t$length\t$genelength\n";
	$long[$n]=$length;
        $genelong[$n]=$genelength;
        $n++;
	}
for ($j=0;$j<$n;$j++)
	{
	$sum=$sum+$long[$j];
        $genesum=$genesum+$genelong[$j];
	}
$average_length=$sum/$n;
$geneaver_length=$genesum/$n;
print "CDS_average\t$average_length\tGENE_average\t$geneaver_length\t$n\n";

sub usage{
print STDERR <<HELP
Usage: perl $0 -i [.GPE] 
        --input|-i:s                     gpe(no bin);
        --bin|-b                         with bin
        --help|-h                        print this message
HELP
} 