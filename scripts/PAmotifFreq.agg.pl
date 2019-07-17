#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use Bio::Perl;
use Bio::DB::Fasta;
my ($motifFile,$refGene,$PAbedfile,$BED,$upstream);
$upstream=50;
GetOptions(
            'm|motif=s'     => \$motifFile,
            'f|fa=s'        => \$refGene,
            'b|bed=s'       => \$PAbedfile,
			'u|up=i'        => \$upstream,
            'h|help'        => sub{usage()}
	  ) || usage();
if(defined $PAbedfile){
	open $BED,"$PAbedfile" or die "Can't open $PAbedfile:$!";
}else{
	$BED=\*STDIN;
}
open MOT,"$motifFile" or die "Can't open $motifFile:$!";
my(%PA,%motif);
my $db = Bio::DB::Fasta->new($refGene);
while(<MOT>){
    chomp;
    $motif{$_}=$_;
}
my $totalPA=0;
while(<$BED>){
    chomp;
	#say STDERR $_;
	$totalPA+=1;
    my @split=split;
    my($chr,$end,$strand)=@split[0,2,5];
    foreach my $motif(keys %motif){
        $motif=~ s/U/T/g;
        if($strand eq '+'){
            my $seq=$db->seq($chr,($end-$upstream),$end+6);
			next if length($seq)<57;
            for(0..$upstream){
                if(substr($seq,$_,6) =~ /$motif/i){
                    $PA{$_-$upstream}{$motif}+=1;
					#say STDERR $_;
                }
            }
        }else{
            my $seq=$db->seq($chr,($end+$upstream),$end-6);
			next if length($seq)<57;
            for(0..$upstream){
                if(substr($seq,$_,6)=~/$motif/i){
                    $PA{$_-$upstream}{$motif}+=1;
					#say STDERR $_;
                }
            }
        }
    }
}
print "pos";
foreach my $motifSeq(sort keys %motif){
	print "\t$motifSeq";
}
print "\n";
foreach my $pos(sort {$a<=>$b} keys %PA){
    print "$pos";
	foreach my $motifSeq(sort keys %{$PA{$pos}}){
		my $freq=sprintf "%.4f",$PA{$pos}{$motifSeq}/$totalPA;
		print "\t$freq";
	}
	print "\n";
}
sub usage{
print <<HELP; 
Usage:	perl $0 -f [reference.fa] -b PA.bed -m motif_seq > result.file
        Statistics the frequency of the input motif in region upstream of PA sites.
Author: Yumei Li, 2018/12/28
	'm|motif=s'     the file storing polyA motif sequence,each at a single line
	'b|bed=s'       the bed file storing PA sites (Can be read from STDIN)
	'u|up=i'        the upstream searching region from PA site[default: 50]
	'f|fa=s'        reference genome, fasta format
	'h|help'        print this help message    
HELP
    exit(-1);
}

