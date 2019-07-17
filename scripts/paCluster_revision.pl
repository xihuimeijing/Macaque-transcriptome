#!/usr/bin/perl -w

=hey
Author: Shijian Sky Zhang
E-mail: zhangsjsky@pku.edu.cn
=cut

use 5.012;
use Getopt::Long;
use File::Basename;
use List::MoreUtils qw(uniq);
use lib dirname $0;
use pm::common;

my ($distance, $winSize, $manner) = (30, 3);
GetOptions(
            'd|distance=i'      => \$distance,
            'm|manner=s'        => \$manner,
            'w|windowSize=i'    => \$winSize,
            'h|help'            => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't open $ARGV[0]: $!";
my (%hash,%geneName); #Revised by Yumei Li, 20181206
while(<IN>){
    chomp;
    my ($chr, $start, $end, $name, $strand, $blockSizes, $blockRelStarts) = (split "\t")[0..3, 5, 10, 11];
    my @blockSizes = split ",", $blockSizes;
    my @blockRelStarts = split ",", $blockRelStarts;
    my ($blockStarts, $blockEnds) = &common::getAbsLoc($start, \@blockSizes, \@blockRelStarts);
    if(exists $hash{$chr}{$strand}){
        push @{$hash{$chr}{$strand}}, [$start, $end, $name];
    }else{
        $hash{$chr}{$strand} = [ [$start, $end, $name] ];
    }
    $geneName{$name}=(split "\t")[-1]; #Revised by Yumei Li, 20181206
}

while(my ($chr, $chrV) = each %hash){
    for my $strand(keys %$chrV){
        my $strandV = $chrV->{$strand};
        my @sortedPAs = $strand eq '+' ? sort {$a->[1]<=>$b->[1]}@$strandV : sort {$a->[0]<=>$b->[0]}@$strandV;
        next if @sortedPAs == 0;
        my ($paRangeStart, $paRangeEnd, @readNames);
        if($strand eq '+'){
            $paRangeEnd = $sortedPAs[0]->[1];
            $paRangeStart = $paRangeEnd - 1;
        }else{
            $paRangeStart = $sortedPAs[0]->[0];
            $paRangeEnd = $paRangeStart + 1;
        }
        @readNames = $sortedPAs[0]->[2];
        my @relPaSites = (1);
        for(my $i = 1; $i <= $#sortedPAs; $i++){
            my $currentPa = $strand eq '+' ? $sortedPAs[$i]->[1] : ($sortedPAs[$i]->[0] + 1);
            if($currentPa - $paRangeEnd <= $distance){
                $paRangeEnd = $currentPa;
                push @readNames, $sortedPAs[$i]->[2];
                push @relPaSites, $currentPa - $paRangeStart;
            }else{
                @relPaSites = reverse map{$paRangeEnd - $paRangeStart - $_ + 1}@relPaSites if $strand eq '-';
                
                # get the mode PA site supporting reads ratio
                my %counts;                
                $counts{$_}++ for(@relPaSites);
                my $modeCount = &modeCountInWindow(\%counts, $winSize);
                
                # get the real PA site
                my @relSite2Count;
                push @relSite2Count, [$_, $counts{$_}] for(keys %counts);
                if($manner eq 'mode'){
                    @relSite2Count = sort{$b->[1]<=>$a->[1] || $b->[0]<=>$a->[0]}@relSite2Count;
                }elsif($manner eq 'downstream'){
                    @relSite2Count = sort{$b->[0]<=>$a->[0]}@relSite2Count;
                }else{ # upstream
                    @relSite2Count = sort{$a->[0]<=>$b->[0]}@relSite2Count;
                }
                my $relPaSite = $relSite2Count[0]->[0];
                my $paSite;
                if($strand eq '+'){
                    $paSite = $paRangeStart + $relPaSite
                }else{
                    @readNames = reverse @readNames;
                    $paSite = $paRangeEnd - $relPaSite + 1;
                }
                
                # get supporting read fraction of each site in cluster
                my @scores = map{0}1..$relPaSites[-1];
                $scores[$_-1]++ for(@relPaSites);
                my $readCount = $#readNames + 1;
                my @freq = map{sprintf "%.2f", $_/$readCount}@scores;
                #Add by Yumei Li
                my @geneNames;
                for(my $j=0;$j<=$#readNames;$j++){
                    push @geneNames,$geneName{$readNames[$j]};
                }
                my @uniqGeneNames=uniq @geneNames;
                #Add by Yumei Li, END 
                say join "\t", ($chr, $paRangeStart, $paRangeEnd,
                                  (join ',', @readNames),
                                  $readCount,
                                  $strand,
                                  $paSite-1, $paSite,
                                  (join ',', @uniqGeneNames),
                                  (join ',', @relPaSites),
                                  scalar @relSite2Count,
                                  (sprintf "%.2f", $modeCount/$readCount),
                                  @freq);
                
                $paRangeStart = $currentPa - 1;
                $paRangeEnd = $currentPa;
                @readNames = ($sortedPAs[$i]->[2]);
                @relPaSites = (1);
            }
        }
        @relPaSites = reverse map{$paRangeEnd - $paRangeStart - $_ + 1}@relPaSites if $strand eq '-';
                
        # get the mode PA site supporting reads ratio
        my %counts;                
        $counts{$_}++ for(@relPaSites);
        my $modeCount = &modeCountInWindow(\%counts, $winSize);
        
        # get the real PA site
        my @relSite2Count;
        push @relSite2Count, [$_, $counts{$_}] for(keys %counts);
        if($manner eq 'mode'){
            @relSite2Count = sort{$b->[1]<=>$a->[1] || $b->[0]<=>$a->[0]}@relSite2Count;
        }else{
            @relSite2Count = sort{$b->[0]<=>$a->[0]}@relSite2Count;
        }
        my $relPaSite = $relSite2Count[0]->[0];
        my $paSite;
        if($strand eq '+'){
            $paSite = $paRangeStart + $relPaSite
        }else{
            @readNames = reverse @readNames;
            $paSite = $paRangeEnd - $relPaSite + 1;
        }
        
        # get supporting read fraction of each site in cluster
        my @scores = map{0}1..$relPaSites[-1];
        $scores[$_-1]++ for(@relPaSites);
        my $readCount = $#readNames + 1;
        my @freq = map{sprintf "%.2f", $_/$readCount}@scores;
        #Add by Yumei Li
        my @geneNames;
        for(my $j=0;$j<=$#readNames;$j++){
                push @geneNames,$geneName{$readNames[$j]};
        }
        my @uniqGeneNames=uniq @geneNames;
        #Add by Yumei Li, END 
        say join "\t", ($chr, $paRangeStart, $paRangeEnd,
                          (join ',', @readNames),
                          $readCount,
                          $strand,
                          $paSite-1, $paSite,
                          (join ',', @uniqGeneNames),
                          (join ',', @relPaSites),
                          scalar @relSite2Count,
                          (sprintf "%.2f", $modeCount/$readCount),
                          @freq);
        
    }
}

sub modeCountInWindow{
    my ($counts, $winSize) = @_;
    my @sites = sort {$a<=>$b}keys %$counts;
    my $currentCount = 0;
    $currentCount += defined $counts->{$sites[0]+$_} ? $counts->{$sites[0]+$_} : 0 for(0..($winSize-1));
    my $modeCount = $currentCount;
    for(my $i = $sites[0]+1; $i <= $sites[-1] - $winSize + 1; $i++){
        $currentCount = $currentCount +
                        (defined $counts->{$i+$winSize-1} ? $counts->{$i+$winSize-1} : 0) -
                        (defined $counts->{$i-1} ? $counts->{$i-1} : 0);
        $modeCount = $currentCount if $currentCount > $modeCount;
    }
    return $modeCount;
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName consensus.bed12 >paCluster.bed8+
Revision: Revised by Yumei Li on 20181206. Change input file to abtain gene information.
    If unambi.bed12+ (generated by readGroup.pl with the last column recording gene name) isn't specified, input from STDIN
    Supporting reads and their number are put into 4th and 5th columns, respectively.
    PA site start and end are put into thickStart and thickEnd columns, respectively.
    Corresponding gene name is put into 9th column.
    Relative positions (comma-separated) of cleavage sites in a cluster and their
        total number are put into 10th and 11th columns, respectively.
    The fraction of supporting reads at the mode cleavage region (in a window, see
        --windowSize) is put into 12th and those at each cleavage sites are appended as additional columns.
Options:
    -d --distance   INT The max distance between two adjacent PA to put them into the same cluster[30]
    -m --manner     STR Manner to pick up PA in cluster. It can be 'mode', 'downstream' or 'upstream'.
    -w --windowSize INT The window size to calculate mode PA site[3].
    -h --help           Print this help information
HELP
    exit(-1);
}