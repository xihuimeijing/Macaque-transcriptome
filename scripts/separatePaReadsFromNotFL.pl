#!/usr/bin/perl -w

=hey
Author: Shijian Sky Zhang
E-mail: zhangsjsky@pku.edu.cn
=cut

use 5.010;
use strict;
use Getopt::Long;
use File::Basename;

my ($polyASize, $windowSize, $fraction) = (20, 50, 0.65);
sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT.fq >OUTPUT.fq
    If INPUT isn't specified, input from STDIN
Options:
    -a --polyA      INT     The minimal polyA length of preliminary trimming[20]
    -w --window     INT     Sliding window size[$polyASize]
    -f --fraction   DOU     Minimal fraction of A base to continue trimming[$fraction]
    -h --help               Print this help information
HELP
    exit(-1);
}

GetOptions(
            'a|polyA=i'     => \$polyASize,
            'w|window=i'    => \$windowSize,
            'f|fraction=s'  => \$fraction,
            'h|help'        => sub{usage()}
        ) || usage();

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't open $ARGV[0]: $!";

while(<IN>){
    chomp;
    my $name = $_;
    chomp(my $seq = <IN>);
    <IN>;
    chomp(my $qual = <IN>);
    
    my $originalSeqLength = length $seq;
    ($seq, $qual) = &trim($seq, $qual);
    my $currentSeqLength = length $seq;
    if($currentSeqLength < $originalSeqLength){
        say $name, "\n$seq\n+\n", $qual if $currentSeqLength > 0;
        say STDERR join "\t", ($name, $originalSeqLength, $currentSeqLength);
    }else{ # no trimming
        $seq = &reverseComplement($seq);
        $qual = join '', reverse (split '', $qual);
        ($seq, $qual) = &trim($seq, $qual);
        $currentSeqLength = length $seq;
        if($currentSeqLength < $originalSeqLength){
            say $name, "\n$seq\n+\n", $qual if $currentSeqLength > 0;
        }
        say STDERR join "\t", ($name, $originalSeqLength, $currentSeqLength);
    }
}

sub trim(){
    my ($seq, $qual) =@_;
    
    # preliminary trimming
    my $preSeqLength = length $seq;
    if($seq =~ s/(A{$polyASize,})$//i){
        $qual = substr $qual, 0, ((length $qual) - (length $1));
    }
    
    my $postSeqLength = length $seq;
    if($postSeqLength >= $windowSize){
        # trim with sliding windows
        my $i = $postSeqLength - $windowSize;
        for(; $i >= 0; $i--){
            last if &getAFraction(substr $seq, $i, $windowSize) < $fraction;
        }
        my $siteTrimTo = $i + $windowSize;
        $seq = substr $seq, 0, $siteTrimTo;
        $qual = substr $qual, 0, $siteTrimTo;
        
        # final trimming
        if(length $seq != $postSeqLength){ # has sliding-window trimming
            if($seq =~ s/(A+)$//i){
                $qual = substr $qual, 0, ((length $qual) - (length $1));
            }
        }
    }
    return ($seq, $qual);
}

sub getAFraction(){
    my ($seq) = @_;
    my $seqLen = length $seq;
    my @As = ($seq =~ /(A)/gi);
    @As / $seqLen;
}

sub reverseComplement(){
    my ($seq) = @_;
    my $revSeq = join '', reverse (split '', $seq);
    $revSeq =~ tr/ATCG/TAGC/;
    $revSeq;
}

