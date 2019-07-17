#!/usr/bin/env perl

use warnings;
use strict;

use Getopt::Long;

use vars qw ($bin $reference $s $n);
use vars qw ($reference );
GetOptions (
        "bin|b"               =>      \$bin,
	"reference|r:s"       =>      \$reference,   
	"name|n"              =>      \$n,
	"strand|s"            =>      \$s,
	"help|h"              =>      sub{&usage;exit(0);}
);

use lib "/rd1/user/chenjy/perl_lib/lib/perl5/";
use Bio::DB::Fasta;
my $ref =$reference;
my $db = Bio::DB::Fasta->new($ref);
my ($line,@str,$strand,$chr,$genestart,$geneend,@exonstart,@exonend,$exonnum,@gene,$gene_id,$trans_id,$seq,$seqN,$num,@base);


while($line=<>)
{
	chomp $line;
	@str=split/\t/,$line;
	shift @str if defined $bin;
	$strand=$str[2];
	$chr=$str[1];
	@exonstart=split/\,/,$str[8];
	@exonend=split/\,/,$str[9];
	$genestart=$str[3];
	$geneend=$str[4];
	$gene_id=$str[11];
	$trans_id=$str[0];
	$exonnum=$str[7];
	
	for (my $i=0;$i<$exonnum;$i++)
		{
		$seq=$seq.$db->seq($chr,$exonstart[$i]+1,$exonend[$i]);
		}
	if (defined $s)
	{
	if ($strand=~/\-/)
		{
		$seqN=reverse $seq;
		$seqN=~tr/ATCGatcg/TAGCtagc/;
		$seq=$seqN;
		}
	}
	@base=split/|/,$seq;
	$num=@base;
	if (defined $n)
	{
	print "\>$gene_id|$trans_id|$chr:$genestart-$geneend|$strand|$num\n";
	}
	print "$seq\n";
	$seq=undef;
#print "@str\n$strand\n$chr\n$start\n$end\n@names\n@gene\n@trans\n";
	}

sub usage{
print STDERR <<HELP
version 20160801 #!perl->#!/usr/bin/env perl
Usage: perl $0 -b -h ***.gpe > ***.fa     
        --bin|-b                    the input with bin column
	--reference|-r              the reference
        --name|n                    print the first line '>name', or only sequence will be printed.
	--strand|-s                 strand-specific, or the strand is the same as reference.
	--help|-h                   print this message
HELP
}
