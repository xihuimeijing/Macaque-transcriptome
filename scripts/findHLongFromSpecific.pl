#!/usr/bin/env perl
use 5.010;
use warnings;
use Getopt::Long;
use File::Basename;
my ($inFile,$RPAfile,$short,$ignoreCase);
GetOptions(
            'i=s'     => \$inFile,
            'm=s'     => \$RPAfile,
	    'g'       => \$ignoreCase,
            's'       => \$short,
            'h|help'  => sub{usage()}
        ) || usage();
my $IN;
if(! defined $inFile){
    $IN=\*STDIN;
}else{
    open $IN,"$inFile" or die "Can't open file $inFile:$!";
}
open RPA,"$RPAfile" or die "Can't open file $RPAfile:$!";
my (%R_PAs,%RgeneInfor);
while(<RPA>){
    chomp;
    my @split=split;
    my $gene=(split /:/,$split[3])[0];
    my $locus=$split[2];
    if(defined $ignoreCase){
    	$gene=uc($gene);
    }
    if(! exists $R_PAs{$gene}{$locus}){
        $R_PAs{$gene}{$locus}=1;
    }else{
        $R_PAs{$gene}{$locus}+=1;
    }
    $RgeneInfor{$gene}{"chr"}=$split[0];
    $RgeneInfor{$gene}{"strand"}=$split[5];
}
while(<$IN>){
    chomp;
    my @split=split;
    my ($gene,$locus)=split /:/,$split[3];
    if(defined $ignoreCase){
            $gene=uc($gene);
    }
    my $H2R_pos=$split[2];
    if(exists $R_PAs{$gene}){
        if($RgeneInfor{$gene}{"chr"} eq $split[0] && $RgeneInfor{$gene}{"strand"} eq $split[5]){
            my ($closestPA,$closestDis,$closestPA_s,$closestDis_s);
            my $tag=0;
            my $tag_short=0;
            if($split[5] eq "+"){
                foreach my $Rpos(keys %{$R_PAs{$gene}}){
                        if($H2R_pos-$Rpos>30){
                                if(! defined $closestPA){
                                    $closestPA=$Rpos;
                                    $closestDis=$H2R_pos-$Rpos;
                                }else{
                                    if($H2R_pos-$Rpos<$closestDis){
                                        $closestPA=$Rpos;
                                        $closestDis=$H2R_pos-$Rpos;
                                    }
                                }
                        }else{
                            $tag+=1;
                        }
                        if(defined $short){
                            if($H2R_pos-$Rpos < -30){
                                if(! defined $closestPA_s){
                                    $closestPA_s=$Rpos;
                                    $closestDis_s=$H2R_pos-$Rpos;
                                }else{
                                    if($H2R_pos-$Rpos>$closestDis_s){
                                        $closestPA_s=$Rpos;
                                        $closestDis_s=$H2R_pos-$Rpos;
                                    }
                                }
                            }else{
                                $tag_short+=1;
                            }
                        }
                }
            }else{
                foreach my $Rpos(keys %{$R_PAs{$gene}}){
                    if($H2R_pos-$Rpos < -30){
                            if(! defined $closestPA){
                                $closestPA=$Rpos;
                                $closestDis=$H2R_pos-$Rpos;
                            }else{
                                if($H2R_pos-$Rpos>$closestDis){
                                    $closestPA=$Rpos;
                                    $closestDis=$H2R_pos-$Rpos;
                                }
                            }
                    }else{
                        $tag+=1;
                    }
                    if(defined $short){
                        if($H2R_pos-$Rpos > 30){
                            if(! defined $closestPA_s){
                                $closestPA_s=$Rpos;
                                $closestDis_s=$H2R_pos-$Rpos;
                            }else{
                                if($H2R_pos-$Rpos<$closestDis_s){
                                    $closestPA_s=$Rpos;
                                    $closestDis_s=$H2R_pos-$Rpos;
                                }
                            }
                        }else{
                            $tag_short+=1;
                        }
                    }
                }
            }
            if($tag == 0){
                say join "\t",$_,$closestPA-1,$closestPA;
            }
            if(defined $short){
                if($tag_short == 0){
                    say STDERR join "\t",$_,$closestPA_s-1,$closestPA_s;
                }
            }
        }
    }
}
sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName -i H.*SpecificPA.Rpos.tsv -m R.*.PA.onGene.bed6 >output.long.bed6+ 2>output.short.bed6+
Author: Yumei Li, 20181216
Revision: Yumei Li, 20181218, Add reporting human specific short PA events.
Revision: Yumei Li, 20190103, Revise to consider chrom and strand information when comparing with macaque PA sites.
Revision: Yumei Li, 20190105, Add parameter to ignore cases when comparing gene names.
Output: With two columns recording the closest upstream/downstream macaque PA site appending to BED6 human specific PA files.
Options:
    -i    FILE      Human specific PA in rheMac8 coordinates in bed6 format (Can be read from STDIN).
    -m    FILE      All macaque PA sites in BED6 format (eg, R.*.PA.onGene.bed6).
    -g    LOGIC     Ignore case when comparing gene names between species.
    -s    LOGIC     Report human short PA events.
    -h --help   Print this help information
HELP
    exit(-1);
}
