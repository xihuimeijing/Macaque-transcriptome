#!/usr/bin/perl -w

=hey
Author: Shijian Sky Zhang
E-mail: zhangsjsky@pku.edu.cn
=cut

use 5.012;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::common;

sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT.bed6 >OUTPUT.bed6+ 2>discarded.bed6+
    If INPUT isn't specified, input from STDIN
Options:
    -b|--twoBitFile     FILE    The 2.bit file
    -t|--paTail         FILE    File with read name and PA tail length
    -a --polyA          INT     The minimal polyA length of preliminary trimming[20]
    -w|--winSize        INT     The window size to scan A base[50]
    -f|--fraction       DOU     The minimal fraction of A base to discard a read[0.5]
    -r|--report         FILE    (Optional) The file to output the statistics
    -h --help                   Print this help information
HELP
    exit(-1);
}

my ($twoBitFile, $paTailFile, $chrSize, $reportFile);
my ($polyASize, $winSize, $fraction) = (20, 50, 0.5);
GetOptions(
            'b|twoBit=s'        => \$twoBitFile,
            't|paTail=s'        => \$paTailFile,
            'a|polyA=i'         => \$polyASize,
            'w|winSize=i'       => \$winSize,
            'f|fraction=f'      => \$fraction,
            'r|report=s'        => \$reportFile,
            'h|help'            => sub{usage()}
        ) || usage();
die "Please specif the 2bit file\n" unless defined $twoBitFile;

$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't read $ARGV[0]: $!";

my %chrSizes;
open CHR, "twoBitInfo $twoBitFile stdout | " or die "Can't read $twoBitFile";
while(<CHR>){
    my ($chr, $size) = split "\t";
    $chrSizes{$chr} = $size;
}

my %paTail;
open TAIL, "$paTailFile" or die "Can't read $paTailFile";
while(<TAIL>){
    chomp;
    my ($name, $length) = split "\t";
    $paTail{$name} = $length;
}

if(defined $reportFile){
    open REPORT, ">$reportFile" or die "Can't write $reportFile";
    say REPORT join "\t", ("#chr", "start", "end", "name", "strand", "Preliminary", "Sliding window", "Final", "On Read");
}

while(<IN>){
    chomp;
    my ($chr, $start, $end, $name, $strand) = (split "\t")[0..3, 5];
    my $chrSize = $chrSizes{$chr};
    my ($genomicAbaseLen, $winStart, $winEnd);
    my ($preliminary, $sliding, $final) = (0);
    my $longEnoughLength = 10_000; # must be larger than $polyASize and $winSize
    
    if($strand eq '+'){
        $longEnoughLength = $chrSize - $end if $end + $longEnoughLength > $chrSize;
        $winStart = $end;
        $winEnd = $winStart + $longEnoughLength;
        my $seq = '';
        if($winStart != $chrSize){
            open FA, "twoBitToFa -seq=$chr -start=$winStart -end=$winEnd $twoBitFile stdout |";
            my @seq = <FA>;
            $seq = join '', map{chomp; $_}@seq[1..$#seq];
        }
        # preliminary trimming
        my $aCount = 0;
        $aCount = length $1 if $seq =~ /(^A+)/i;
        if($aCount >= $polyASize){
            $winStart += $aCount;
            $preliminary = $winStart - $end;
            $seq =~ s/(^A+)//i;
        }
        # sliding window
        for(my $i = 0; $i <= (length $seq) - $winSize; $i++){
            my @As = ((substr $seq, $i, $winSize) =~ /(A)/ig);
            if(@As/$winSize < $fraction){
                $winStart += $i;
                $seq = substr $seq, $i;
                last;
            }
        }
        $sliding = $winStart - $end;
        # final trimming
        $winStart += length $1 if $seq =~ s/(^A+)//i;
        $genomicAbaseLen = $winStart - $end;
    }else{
        $longEnoughLength = $start if $start - $longEnoughLength < 0;
        $winEnd = $start;
        $winStart = $winEnd - $longEnoughLength;
        my $seq = '';
        if($winEnd != 0){
            open FA, "twoBitToFa -seq=$chr -start=$winStart -end=$winEnd $twoBitFile stdout |";
            my @seq = <FA>;
            $seq = join '', map{chomp; $_}@seq[1..$#seq];
        }
        # preliminary trimming
        my $tCount = 0;
        $tCount = length $1 if $seq =~ /(T+$)/i;
        if($tCount >= $polyASize){
            $winEnd -= $tCount;
            $preliminary = $start - $winEnd;
            $seq =~ s/(T+$)//i;
        }
        # sliding window
        for(my $i = length($seq); $i - $winSize >= 0; $i--){
            my @Ts = ((substr $seq, $i - $winSize, $winSize) =~ /(T)/ig);
            if(@Ts/$winSize < $fraction){
                $winEnd -= (length $seq) - $i;
                $seq = substr $seq, 0, $i;
                last;
            }
        }
        $sliding = $start - $winEnd;
        # final trimming
        $winEnd -= length $1 if $seq =~ s/(T+$)//i;
        $genomicAbaseLen = $start - $winEnd;
    }
    if(exists $paTail{$name}){
        if($paTail{$name} - $genomicAbaseLen > 20){
            say;
        }else{
            say STDERR;
        }
        say REPORT join "\t", ($chr, $start, $end, $name, $strand, $preliminary, $sliding, $genomicAbaseLen, $paTail{$name}) if defined $reportFile;
    }else{
        say STDERR;
        say REPORT join "\t", ($chr, $start, $end, $name, $strand, $preliminary, $sliding, $genomicAbaseLen, "NA") if defined $reportFile;
    }
}
