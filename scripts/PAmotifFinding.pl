#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use Bio::Perl;
use Bio::DB::Fasta;
my ($motifFile,$refGene,$PAbedfile,$length,$BED);
$length=100;
GetOptions(
            'f|fa=s'        => \$refGene,
            'b|bed=s'       => \$PAbedfile,
            'l|len=i'		=> \$length,
			'h|help'        => sub{usage()}
	  ) || usage();
$motifFile='/mnt/share/liym/data/PAmotif.12.txt';
if(! defined $PAbedfile){
	$BED=\*STDIN;
}else{
	open $BED,"$PAbedfile" or die "Can't open $PAbedfile:$!";
}
open MOT,"$motifFile" or die "Can't open $motifFile:$!";
my(%PA,@motif);
my $db = Bio::DB::Fasta->new($refGene);
while(<MOT>){
    chomp;
    push @motif,$_;
}
while(<$BED>){
    chomp;
    my @split=split;
    my($chr,$start,$end,$strand)=@split[0,1,2,5];
	my $output=0;
	LINE:for(my $i=0;$i<=$#motif;$i++){
        $motif[$i]=~ s/U/T/g;
        if($strand eq '+'){
            my $seq=$db->seq($chr,($end-$length),$end);
            for(0..($length-5)){
                if(substr($seq,$_,6) =~ /$motif[$i]/i){
                    say join "\t",@split,$motif[$i],$start-$length+$_,$start-$length+$_+6,$_;
					$output+=1;
                    last LINE;
                }
            }
        }else{
            my $seq=$db->seq($chr,($end+$length),$end);
            for(0..($length-5)){
                if(substr($seq,$_,6)=~/$motif[$i]/i){
                    say join "\t",@split,$motif[$i],$start+$_,$start+$_+6,$_;
					$output+=1;
                    last LINE;
                }
            }
        }
    }
	if($output == 0){
		say join "\t",@split,"NA","NA","NA","NA";
	}
}
sub usage{
print <<HELP; 
Usage:	perl $0 -f [reference.fa] -b PA.bed > result.file
        Output the PA motif type and coordinate for each line.
Output: bed6_fields PAmotif motif_start motif_end motifDistanceToEnd
Author: Yumei Li, 20181108
Revision: Yumei Li, 20181213, Revised to output only one motif for each PA site.
		'b|bed=s'       the bed file storing PA sites
		'l|len=i'       the upstream length to search PA motif [default: 100].
		'f|fa=s'        reference genome, fasta format
		'h|help'        print this help message    
HELP
    exit(-1);
}
