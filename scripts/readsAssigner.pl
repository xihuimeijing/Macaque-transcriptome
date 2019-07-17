#!/bin/env perl

=hey
Author: Shijian Sky Zhang
E-mail: zhangsjsky@pku.edu.cn
=cut

use 5.012;
use warnings;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::common;

my ($gpeFile, $bin, $singleLine);
my ($offset, $minConsensusIntronN, $minCovOnRead) = (10, 1, 0.9);
my $useStrand;

GetOptions(
            'g|gpe=s'               => \$gpeFile,
            'b|bin'                 => \$bin,
            'o|offset=i'            => \$offset,
            'i|consensusIntron=i'   => \$minConsensusIntronN,
            'r|minCoverageOnRead=s' => \$minCovOnRead,
            's|singleLine'          => \$singleLine,
            'strand'                => \$useStrand,
            'h|help'                => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't open $ARGV[0]: $!";
my $GPE;
open $GPE, "$gpeFile" or die "Can't open $gpeFile: $!";

my %hash = &common::buildHash2($GPE, $bin);

while(<IN>){
    chomp;
    my ($chr, $start, $end, $strand, $blockSizes, $relStarts) = (split "\t")[0, 1, 2, 5, 10, 11];
    my @relStarts = split ",", $relStarts;
    my @readSizes = split ",", $blockSizes;
    my ($readStarts, $readEnds) = &common::getAbsLoc($start, \@readSizes, \@relStarts);
    my (@matchTrans, @matchCoverageOnRead, @matchCoverageOnTrans, @matchTransExonCount, @matchGenes);
    if(@readSizes == 1){ # single-exon read
        my @inEndExonAndExtension;
        my @strandToCheck = defined $useStrand ? ($strand) : ('+', '-');
        for my $strand (@strandToCheck){
            my $strandV = $hash{$chr}->{$strand};
            for my $trans (@$strandV){
                if($start < $trans->[1] && $end > $trans->[0]){ # intersect
                    my ($transName, $transBlockStarts, $transBlockEnds, $geneName) = (split "\t", $trans->[2])[0, 8, 9, 11];
                    my @transBlockStarts = split ',', $transBlockStarts;
                    my @transBlockEnds = split ',', $transBlockEnds;
                    my $overlapLen = &common::getOverlapLength(\@transBlockStarts, \@transBlockEnds, $readStarts, $readEnds);
                    $overlapLen = 0 unless defined $overlapLen;
                    my $coverageOnRead = $overlapLen / &common::getExonsLength($readStarts, $readEnds);
                    my $coverageOnTrans = $overlapLen / &common::getExonsLength(\@transBlockStarts, \@transBlockEnds);
                    push @matchTrans, $transName;
                    push @matchCoverageOnRead, $coverageOnRead;
                    push @matchCoverageOnTrans, $coverageOnTrans;
                    push @matchTransExonCount, $#transBlockStarts + 1;
                    push @matchGenes, $geneName;
                    if($transBlockStarts[0] < $end && $end <= $transBlockEnds[0] ||
                       $transBlockStarts[-1] <= $start && $start < $transBlockEnds[-1]){
                        push @inEndExonAndExtension, 1;
                    }else{
                        push @inEndExonAndExtension, 0;
                    }
                }
                last if $trans->[0] > $end;
            }
        }
        if(@matchTrans == 0){
            say join "\t", ($_, 'IG', 'NA', 'NA', 'NA', 'NA', 'NA');
        }else{
            if(defined $singleLine){
                my $type = 'I';
                my (@realMatchTrans, @realMatchCoverageOnRead, @realMatchCoverageOnTrans, @realMatchGenes);
                for(my $i = 0; $i <= $#matchTrans; $i++){
                    next if $matchTransExonCount[$i] > 1 && $inEndExonAndExtension[$i] == 0 && $matchCoverageOnRead[$i] < $minCovOnRead;
                    $type = 'E';
                    push @realMatchTrans, $matchTrans[$i];
                    push @realMatchCoverageOnRead, $matchCoverageOnRead[$i];
                    push @realMatchCoverageOnTrans, $matchCoverageOnTrans[$i];
                    push @realMatchGenes, $matchGenes[$i];
                }
                if($type eq 'E'){
                    say join "\t", ($_, $type,
                                    (join ',', @realMatchTrans),
                                    (join ',', @realMatchGenes),
                                    'NA',
                                    (join ',', @realMatchCoverageOnRead),
                                    (join ',', @realMatchCoverageOnTrans));
                }else{
                    say join "\t", ($_, $type,
                                    (join ',', @matchTrans),
                                    (join ',', @matchGenes),
                                    'NA',
                                    (join ',', @matchCoverageOnRead),
                                    (join ',', @matchCoverageOnTrans));
                }
            }else{
                for(my $i = 0; $i <= $#matchTrans; $i++){
                    my $type = $matchTransExonCount[$i] > 1 && $inEndExonAndExtension[$i] == 0 && $matchCoverageOnRead[$i] < $minCovOnRead ? 'I' : 'E';
                    say join "\t", ($_, $type, $matchTrans[$i], $matchGenes[$i], 'NA', $matchCoverageOnRead[$i], $matchCoverageOnTrans[$i]);
                }
            }
        }
    }else{ # mult-exon read
        my (@consensusIntronN, @juncChainFlank);
        for my $trans (@{$hash{$chr}{$strand}}){
            if($start < $trans->[1] && $end > $trans->[0]){ # intersect
                my ($transName, $transBlockStarts, $transBlockEnds, $geneName) = (split "\t", $trans->[2])[0, 8, 9, 11];
                my @transBlockStarts = split ',', $transBlockStarts;
                my @transBlockEnds = split ',', $transBlockEnds;
                my $overlapLen = &common::getOverlapLength(\@transBlockStarts, \@transBlockEnds, $readStarts, $readEnds);
                $overlapLen = 0 unless defined $overlapLen;
                my $coverageOnRead = $overlapLen / &common::getExonsLength($readStarts, $readEnds);
                my $coverageOnTrans = $overlapLen / &common::getExonsLength(\@transBlockStarts, \@transBlockEnds);
                my $consensusIntronN = &common::getConsensusIntronN(\@transBlockStarts, \@transBlockEnds, $readStarts, $readEnds, $offset);
                my $transJuncChain = join ';', map{"$transBlockEnds[$_]-$transBlockStarts[$_+1]"}0..($#transBlockStarts-1);
                my $readJuncChain = join ';', map{"$readEnds->[$_]-$readStarts->[$_+1]"}0..(@$readStarts-2);
                push @juncChainFlank, ($transJuncChain =~ /$readJuncChain$/ || $transJuncChain =~ /^$readJuncChain/ ? 1 : 0);
                push @matchTrans, $transName;
                push @matchCoverageOnRead, $coverageOnRead;
                push @matchCoverageOnTrans, $coverageOnTrans;
                push @matchGenes, $geneName;
                push @consensusIntronN, $consensusIntronN;
            }
            last if $trans->[0] > $end;
        }
        if(@matchTrans == 0){
            say join "\t", ($_, 'IG', 'NA', 'NA', 'NA', 'NA', 'NA');
        }else{
            my $readIntronN = @readSizes - 1;
            if(defined $singleLine){
                my $type = 'I';
                my (@realMatchTrans, @realConsensusIntroN, @realMatchCoverageOnRead, @realMatchCoverageOnTrans, @realMatchGenes);
                for(my $i = 0; $i <= $#matchTrans; $i++){
                    if($juncChainFlank[$i] == 1 ||
                       ($readIntronN < $minConsensusIntronN && $readIntronN == $consensusIntronN[$i]) ||
                        ($readIntronN >= $minConsensusIntronN && $consensusIntronN[$i] >= $minConsensusIntronN)    ){
                        $type = 'E';
                        push @realMatchTrans, $matchTrans[$i];
                        push @realConsensusIntroN, $consensusIntronN[$i];
                        push @realMatchCoverageOnRead, $matchCoverageOnRead[$i];
                        push @realMatchCoverageOnTrans, $matchCoverageOnTrans[$i];
                        push @realMatchGenes, $matchGenes[$i];
                    }
                }
                if($type eq 'E'){
                    say join "\t", ($_, $type,
                                    (join ',', @realMatchTrans),
                                    (join ',', @realMatchGenes),
                                    (join ',', @realConsensusIntroN),
                                    (join ',', @realMatchCoverageOnRead),
                                    (join ',', @realMatchCoverageOnTrans));
                }else{
                    say join "\t", ($_, $type,
                                    (join ',', @matchTrans),
                                    (join ',', @matchGenes),
                                    (join ',', @consensusIntronN),
                                    (join ',', @matchCoverageOnRead),
                                    (join ',', @matchCoverageOnTrans));
                }
            }else{
                for(my $i = 0; $i <= $#matchTrans; $i++){
                    my $type = 'I';
                    if($juncChainFlank[$i] == 1 ||
                       (@readSizes < $minConsensusIntronN && @readSizes-1 == $consensusIntronN[$i]) ||
                       (@readSizes >= $minConsensusIntronN && $consensusIntronN[$i] >= $minConsensusIntronN)){
                        $type = 'E';
                    }
                    say join "\t", ($_, $type, $matchTrans[$i], $matchGenes[$i], $consensusIntronN[$i], $matchCoverageOnRead[$i], $matchCoverageOnTrans[$i]);
                }
            }
        }
    }
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT.bed12 >OUTPUT.bed12+
    If INPUT isn't specified, input from STDIN
    OUTPUT.tsv contains 5 more additional columns:
        For no -s: type, matched transcript, matched gene, consensus intron number, coverage on read, coverage on transcript.
        For -s: matched transcript, matched gene, consensus intron number, coverage on read, coverage on transcript are comma-separated items.
        type may be 'E' (exonic), 'I' (intronic) or 'IG' (intergenic).
Options:
    -g --gpe                FILE    Gene annotation file in gpe format
    -b --bin                        With bin column in --gpe
    -o --offset             INT     Offset for intron boundaries used for finding consensus intron[10]
    -i --consenesusIntron   INT     The min consensus intron number to assign a read to a transcript[1]
    -r --minCoverageOnRead  DOU     The min coverage on read between the read and a transcript to assign the read to the transcript[0.9]
    -s --singleLine                 Output matchs in one line per read
       --strand                     For single-exon read, take the strand column as its real strand
    -h --help                       Print this help information
HELP
    exit(-1);
}