#!/usr/bin/env perl

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::gpeParser;
use pm::bedParser;

GetOptions(
		'h|help'=>	=> sub{usage()}
                )||usage(); 

$ARGV[0]='-' unless defined $ARGV[0];
open BED,"$ARGV[0]" or die "Can't open $ARGV[0]:$!";

while(<BED>){
    chomp;
    my ($chr, $start, $name, $strand, $cdsStart, $cdsEnd, $sizes, $relStarts) = (split "\t")[0, 1, 3, 5..7, 10, 11];
    my @sizes = split ',', $sizes;
    my @relStarts = split ',', $relStarts;
    my ($blockStarts, $blockEnds) = bedParser::getAbsLoc($start, \@sizes, \@relStarts);
    my $exonFrames = join ',', @{&gpeParser::getExonFrames($cdsStart, $cdsEnd, $strand, (join ',', @$blockStarts), (join ',', @$blockEnds))};
    say join "\t",($name, $chr, $strand, $blockStarts->[0], $blockEnds->[-1], $cdsStart, $cdsEnd, $#sizes+1,
                          (join ',', @$blockStarts), (join ',', @$blockEnds), 0, "gene name", "unk", "unk", $exonFrames);
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName Input.bed >Output.gpe
    If Input.bed not specified, input from STDIN
    Output to STDOUT

    -h --help	Print this help information
HELP
    exit(-1);
}