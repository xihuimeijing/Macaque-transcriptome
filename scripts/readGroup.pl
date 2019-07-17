#!/usr/bin/env perl

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

my ($gpeFile, $bin);
GetOptions(
            'g|gpe=s'   => \$gpeFile,
            'b|bin'     => \$bin,
            'h|help'    => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't read file ($ARGV[0]): $!";
open GPE, $gpeFile or die "Can't read file ($gpeFile): $!";

my %trans2Gene;
while(<GPE>){
    chomp;
    my @fields = split "\t";
    shift @fields if defined $bin;
    my ($trans, $gene) = @fields[0, 11];
    $trans2Gene{$trans} = $gene;
}

my %novelReads;
READ:
while(<IN>){
    chomp;
    my @fields = split "\t";
    my ($chr, $strand) = @fields[0, 5];
    my @transs = split ',', $fields[15];
    if($transs[0] eq 'NA' || $fields[14] ne 'E'){
        push @{$novelReads{"$chr:$strand"}}, \@fields;
    }else{
        my @genes = map{$trans2Gene{$_}}@transs;
        my @uniqGenes = &common::uniqArray(\@genes);
        if(@uniqGenes == 1){
            say join "\t", ($_, join ',', @uniqGenes);
        }else{
            say STDERR join "\t", ($_, join ',', @uniqGenes);
        }
    }
}

while(my ($chrStrand, $chrStrandV) = each %novelReads){
    my @sortedStrandV = sort {$a->[1]<=>$b->[1]}@{$chrStrandV};
    my $clusterEnd = $sortedStrandV[0]->[2];
    my $inc = 1;
    say join "\t", (@{$sortedStrandV[0]}, "$chrStrand:$inc");
    for(my $i = 1; $i <= $#sortedStrandV; $i++){
        my ($start, $end) = @{$sortedStrandV[$i]}[1,2];
        if($start < $clusterEnd){
            $clusterEnd = $end if $end > $clusterEnd;
        }else{
            $inc++;
            $clusterEnd = $end;
        }
        say join "\t", (@{$sortedStrandV[$i]}, "$chrStrand:$inc");
    }
}

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT >OUTPUT
    If INPUT isn't specified, input from STDIN
Options:
    -g --gpe    FILE    The gene annotation file in gpe format
    -b --bin            With bin column in --gpe
    -h --help           Print this help information
HELP
    exit(-1);
}