#!/usr/bin/perl -w

=hey
Author: Shijian Sky Zhang
E-mail: zhangsjsky@pku.edu.cn
=cut

use strict;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::common;

my ($gpeFile, $bin, $juncFile);
my ($minRatio, $maxRatio) = (0.5, 2);
GetOptions(
            'g|gpe=s'       => \$gpeFile,
            'b|bin'         => \$bin,
            'j|junc=s'      => \$juncFile,
            'i|minRatio=s'  => \$minRatio,
            'm|maxRatio=s'  => \$maxRatio,
            'h|help'        => sub{usage()}
        ) || usage();

my ($GPE, $JUNC);
open $GPE, "$gpeFile" or die "Can't open $gpeFile: $!";
open $JUNC, "$juncFile" or die "Can't open $juncFile: $!";
$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't open $ARGV[0]: $!";

my %refHash = &common::buildHash2($GPE, $bin);
my %juncHash = &common::buildBedHash($JUNC);

while(<IN>){
    chomp;
    my @fields = split "\t";
    my ($chr, $start, $end, $strand, $readSizes, $relStarts) = @fields[0..2, 5, 10, 11];
    my @readSizes = split ",", $readSizes;
    if(@readSizes == 1){
        print "$_\n";
        next;
    }
    my @relStarts = split ",", $relStarts;
    my ($readStarts, $readEnds) = &common::getAbsLoc($start, \@readSizes, \@relStarts);
    my @readStartsCopy = @$readStarts;
    my @readEndsCopy = @$readEnds;
    my (@trans, @juncs);
    for my $trans (@{$refHash{$chr}{$strand}}){
        if($start < $trans->[1] && $end > $trans->[0]){ # intersect
            push @trans, $trans->[2];
        }
        last if $trans->[0] > $end;
    }
    for my $junc(@{$juncHash{$chr}{$strand}}){
        if($start < $junc->[1] && $end > $junc->[0]){ # intersect
            push @juncs, [ $junc->[2], $junc->[3] ];
        }
        last if $junc->[0] > $end;
    }
    if(@trans > 0 || @juncs > 0){
        for(my $i = 1; $i < $#readSizes; $i++){
            my ($exonStart, $exonEnd) = ($readStarts->[$i], $readEnds->[$i]);
            my ($hasOverlapTransExon, $hasOverlapJuncExon) = (0, 0);
            my (@intronStarts, @intronEnds);
    TRANS:  for my $trans (@trans){
                my ($transStarts, $transEnds) = (split "\t", $trans)[8, 9];
                my @transStarts = split ",", $transStarts;
                my @transEnds = split ",", $transEnds;
                for(my $j = 0; $j <= $#transStarts; $j++){
                    if($transStarts[$j] > $exonEnd){
                        if(defined $transStarts[$j-1] &&
                           $readStarts->[$i-1] < $transEnds[$j-1] && $readEnds->[$i-1] > $transStarts[$j-1] &&
                           $readStarts->[$i+1] < $transEnds[$j] && $readEnds->[$i+1] > $transStarts[$j]){
                            push @intronStarts, $transEnds[$j-1];
                            push @intronEnds, $transStarts[$j];
                        }
                        last;
                    }
                    if($exonStart < $transEnds[$j] && $exonEnd > $transStarts[$j]){
                        $hasOverlapTransExon = 1;
                        last TRANS;
                    }
                }
            }
            if($hasOverlapTransExon == 0){
        JUNC:   for my $junc (@juncs){
                    my ($juncExonStarts, $juncExonEnds) = @$junc;
                    for(my $j = 0; $j <= @$juncExonStarts -1; $j++){
                        if($juncExonStarts->[$j] > $exonEnd){
                            if(defined $juncExonStarts->[$j-1] &&
                               $readStarts->[$i-1] < $juncExonEnds->[$j-1] && $readEnds->[$i-1] > $juncExonStarts->[$j-1] &&
                               $readStarts->[$i+1] < $juncExonEnds->[$j] && $readEnds->[$i+1] > $juncExonStarts->[$j]){
                                push @intronStarts, $juncExonEnds->[$j-1];
                                push @intronEnds, $juncExonStarts->[$j];
                            }
                            last;
                        }
                        if($exonStart < $juncExonEnds->[$j] && $exonEnd > $juncExonStarts->[$j]){
                            $hasOverlapJuncExon = 1;
                            last JUNC;
                        }
                    }
                }
                if($hasOverlapJuncExon == 0){
                    if(&common::uniqArray(\@intronStarts) == 1 && &common::uniqArray(\@intronEnds) == 1){
                        my $diffLen = ($intronStarts[0] - $readEnds->[$i-1]) + ($readStarts->[$i+1] - $intronEnds[0]);
                        my $diffRatio =$diffLen / ($exonEnd - $exonStart);
                        if($minRatio <= $diffRatio && $diffRatio <= $maxRatio){
                            ($readStartsCopy[$i], $readEndsCopy[$i]) = (undef, undef);
                            $readStartsCopy[$i+1] = $intronEnds[0];
                            $readEndsCopy[$i-1] = $intronStarts[0];
                        }
                    }
                }
            }
            
        }
        @readStartsCopy = &common::removeUndefElement(\@readStartsCopy);
        @readEndsCopy = &common::removeUndefElement(\@readEndsCopy);
        @fields[9..11] = ( scalar @readStartsCopy,
                           (join ',', &common::getSizes(\@readStartsCopy, \@readEndsCopy)),
                           (join ',', &common::getRelStarts(\@readStartsCopy))
                         );
        print join "\t", @fields;
        print "\n";
    }else{
        print $_;
        print "\n";
    }
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT >OUTPUT
    If INPUT isn't specified, input from STDIN
Options:
    -g --gpe        FILE    Gene annotaion file in gpe format
    -b --bin                With bin column
    -j --junc       FILE    The junctions.bed file generated from RNA-seq of the PacBio tissue
    -i --minRatio   DOU     The min ratio of a exon length and its flanking mismaped sequence. Exon with ratio between --minRatio and --maxRatio will be removed and the flanking sequence will be refined on the basis of annotation from --gpe[0.5]
    -a --maxRatio   DOU     The max ratio (refer to --minRatio)[2]
    -h --help               Print this help information
HELP
    exit(-1);
}