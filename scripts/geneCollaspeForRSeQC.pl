#!/usr/bin/env perl
use 5.010;
use warnings;
use Getopt::Long;
use File::Basename;
use List::MoreUtils qw(uniq);
my ($inFile);
GetOptions(
            'i=s'     => \$inFile,
            'h|help'    => sub{usage()}
        ) || usage();
my ($IN,%exon,%gene);
if(defined $inFile){
    open $IN,"$inFile" or die "Can't open file $inFile:$!";
}else{
    $IN=\*STDIN;
}
while(<$IN>){
    chomp;
    my @split=split;
    my ($chr,$strand,$txStart,$txEnd,$starts,$ends,$geneName)=@split[1,2,5,6,8,9,11];
    #Remove non-coding transcripts
    next if ($txStart == $txEnd);
    my @starts=split /,/,$starts;
    my @ends=split /,/,$ends;
    for(my $i=0;$i<=$#starts;$i++){
        $gene{$chr}{$strand}{$geneName}{$starts[$i]}{$ends[$i]}="";
        if(! exists $exon{$chr}{$strand}{$starts[$i]}{$ends[$i]}){
            $exon{$chr}{$strand}{$starts[$i]}{$ends[$i]}=$geneName;
        }else{
            my @genes=split /;/,$exon{$chr}{$strand}{$starts[$i]}{$ends[$i]};
            push @genes,$geneName;
            my @uniqGenes=uniq @genes;
            $exon{$chr}{$strand}{$starts[$i]}{$ends[$i]}=join ";",@uniqGenes;
        }
    }
}
#Record regions overlapping with more than one gene
my (%outGPE,%overlapRegions);
my($recordChr,$recordStrand,$recordStart,$recordEnd,$recordGene);
foreach my $chr(keys %exon){
    foreach my $strand(keys %{$exon{$chr}}){
        foreach my $start(sort{$a<=>$b} keys %{$exon{$chr}{$strand}}){
            foreach my $end(sort{$a<=>$b} keys %{$exon{$chr}{$strand}{$start}}){
                if($exon{$chr}{$strand}{$start}{$end} =~ /,/){
                    $overlapRegions{$chr}{$strand}{$start}{$end}="";     
                }else{
                    if(!defined $recordChr){
                        ($recordChr,$recordStrand,$recordStart,$recordEnd)=($chr,$strand,$start,$end);
                        $recordGene=$exon{$chr}{$strand}{$start}{$end};
                    }else{
                        if($chr eq $recordChr && $strand eq $recordStart){
                            if($start>=$recordStart && $start<=$recordEnd && $exon{$chr}{$strand}{$start}{$end} ne $recordGene){#Overlap with different genes
                                    my $newEnd=$recordEnd>$end ? $recordEnd : $end;
                                    $overlapRegions{$chr}{$strand}{$start}{$newEnd}="";
                            }
                        }
                    }
                }
            }
        }
    }
}
#Merge intervals for each gene
foreach my $chr(keys %gene){
    foreach my $strand(keys %{$gene{$chr}}){
        foreach my $geneName(sort keys %{$gene{$chr}{$strand}}){
            undef $recordChr;
            undef $recordStrand;
            undef $recordStart;
            undef $recordEnd;
            undef $recordGene;
            foreach my $start(sort {$a<=>$b} keys %{$gene{$chr}{$strand}{$geneName}}){
                foreach my $end(sort {$a<=>$b} keys %{$gene{$chr}{$strand}{$geneName}{$start}}){
                    if(!defined $recordChr){
                        ($recordChr,$recordStrand,$recordStart,$recordEnd)=($chr,$strand,$start,$end);
                    }else{
                        if($start>=$recordStart && $start<=$recordEnd){
                            my $newEnd=$recordEnd>$end ? $recordEnd : $end;
                            delete $gene{$chr}{$strand}{$geneName}{$recordStart}{$recordEnd}; #Remove regions that overlapped 
                            delete $gene{$chr}{$strand}{$geneName}{$start}{$end};
                            $recordEnd=$newEnd;
                        }else{
                            if(! exists $gene{$chr}{$strand}{$geneName}{$recordStart}{$recordEnd}){
                                $gene{$chr}{$strand}{$geneName}{$recordStart}{$recordEnd}="";
                            }
                        }      
                    }   
                }
                if(! exists $gene{$chr}{$strand}{$geneName}{$recordStart}{$recordEnd}){
                    $gene{$chr}{$strand}{$geneName}{$recordStart}{$recordEnd}="";
                }
            }
        }
        
    }
}
#Remove regions that overlapped with more than one gene
foreach my $chr(keys %gene){
    foreach my $strand(keys %{$gene{$chr}}){
        foreach my $geneName(sort keys %{$gene{$chr}{$strand}}){
            my ($exonS,$exonE);
            foreach my $start(sort {$a<=>$b} keys %{$gene{$chr}{$strand}{$geneName}}){
                my $end=keys %{$gene{$chr}{$strand}{$geneName}{$start}};
                #remove
                foreach my $OregionS(sort {$a<=>$b} keys %{$overlapRegions{$chr}{$strand}}){
                    last if $OregionS>$end;
                    foreach my $OregionE(sort {$a<=>$b} keys %{$overlapRegions{$chr}{$strand}{$OregionS}}){
                        if($OregionE>=$start && $OregionE<=$end){
                            my($interval1_S,$interval1_E)=($start,$OregionS);
                            my($interval2_S,$interval2_E)=($OregionE,$end);
                            if($interval1_E>$interval1_S){
                                if(! defined $exonS){
                                    $exonS=$interval1_S;
                                    $exonE=$interval1_E;
                                }else{
                                    $exonS=$exonS.",".$interval1_S;
                                    $exonE=$exonE.",".$interval1_E;
                                }
                            }
                             if($interval2_E>$interval2_S){
                                if(! defined $exonS){
                                    $exonS=$interval2_S;
                                    $exonE=$interval2_E;
                                }else{
                                    $exonS=$exonS.",".$interval2_S;
                                    $exonE=$exonE.",".$interval2_E;
                                }
                            }
                        }
                    }
                }
            }
            if(defined $exonS){
                my @outStarts=split /,/,$exonS;
                my @outEnds=split /,/,$exonE;
                my @outshift=(-1) x ($#outEnds+1);
                say join "\t",$geneName,$chr,$strand,$outStarts[0],$outEnds[$#outEnds],$outStarts[0],$outEnds[$#outEnds],$#outEnds+1,$exonS,$exonE,0,$geneName,"none","none",(join ",",@outshift);
            }
        }
    }
}
sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName -i <IN.gpe> 
Description: This script collapsed all transcripts to a single transcript model for each gene based on procedures using by GTEx (https://gtexportal.org/home/documentationPage#staticTextAnalysisMethods) 
    (1) Exons associated with non-coding transcripts were excluded.
    (2) Exon intervals overlapping within a gene were merged.
    (3) The intersections of exon intervals overlapping between genes were excluded.
    (4) The remaining exon intervals were mapped to their respective gene identifier in GPE format.
Author: Yumei Li, 2018/12/20
Options:
    -i    FILE  Gene structure in GPE format.
    -h --help   Print this help information
HELP
    exit(-1);
}