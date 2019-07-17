#!/usr/bin/env perl
use 5.010;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use lib dirname $0;
use pm::fqParser;

my ($minLength, $repeat, $baseFrac, $qualCutoff, $qualFrac, $minAveQual, $tail, $discardFile);
my ($base, $platform) = ('N', 'Sanger');

sub usage{
   my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName -f 0.1 -q 0.5 -c 5 input.fq[.gz] >output.fq 2>report.log
      If input.fq[.gz] isn't specified, input from STDIN
      Statistic report is outputted to STDERR

   Sequence filter
      -f|--baseFrac     FLOAT       Discard reads with fraction of ambiguous bases > FLOAT. 0.1 is recommended
      -b|--base         CHAR        Ambiguous base[N]
      -l|--minLength    INT         Discard reads with length < INT
      -r|--repeat       INT         Discard reads with poly A, T, C or G in length >= INT
      
   Quality filter
      -q|--qualFrac     FLOAT       Discard reads with fraction of low-quality bases > FLOAT. 0.5 is recommended
      -c|--qualCutoff   INT         Low-quality base when its quality < INT. 5 is recommended
      -a|--minAveQual   FLOAT       Discard reads with average quality < FLOAT
      -t|--tail         INT1,INT2   Read is discarded when any base quality of tail INT1 ones < INT2
                                    E.g.: '2,5' means discarding read when any base quality of tail 2 ones is less than 5
      -p|--platform     CHAR|STR    The quality scoring platform. It can be
                                    [S|Sanger] (Phred+33, 0~40),
                                    X|Solexa (Solexa+64, -5~40),
                                    I|Illumina1.3+ (Phred+64, 0~40),
                                    J|Illumina1.5+ (Phred+64, 3~40),
                                    L|Illumina1.8+ (Phred+33, 0~41)
   Other
      -d|--discard      FILE        Output discarded reads to FILE
      -h|--help                     This help information
HELP
    exit -1;
}

GetOptions(
                        'f|baseFrac=s'    => \$baseFrac,
                        'b|base=s'        => \$base,
                        'l|minLength=i'   => \$minLength,
                        'r|repeat=i'      => \$repeat,
                        'q|qualFrac=s'    => \$qualFrac,
                        'c|qualCutoff=s'  => \$qualCutoff,
                        'a|aveQual=s'     => \$minAveQual,
                        't|tail=s'        => \$tail,
                        'p|platform=s'    => \$platform,
                        'd|discard=s'     => \$discardFile,
                        'h|help'          => sub{usage()}
) || usage();

if(defined $ARGV[0]){
    if(`file -L $ARGV[0]` =~ /gzip/){
        open IN,"gzip -dc $ARGV[0]|" or die "Can't open $ARGV[0]: $!";
    }else{
        open IN,"$ARGV[0]" or die "Can't open $ARGV[0]: $!";
    }
}else{
    open IN, "-";
}
open DIS, ">$discardFile" or die "Can't open $discardFile: $!" if defined $discardFile;
my $readName;
my ($passedN, $badTailN, $totalN) = (0, 0, 0);

while($readName=<IN>){
    chomp(my $seq=<IN>);
    my $plus=<IN>;
    chomp(my $qualSeq=<IN>);
    $totalN++;
    my $toOutput = 1;
    if(defined $minLength){
        $toOutput = 0 if(length $seq < $minLength);
    }
    if($toOutput ==1 && defined $repeat){
        my ($repeatAs, $repeatTs, $repeatCs, $repeatGs) = ('A' x $repeat, 'T' x $repeat, 'C' x $repeat, 'G' x $repeat);
        $toOutput = 0 if $seq =~ /($repeatAs)|($repeatTs)|($repeatCs)|($repeatGs)/;
    }
    if($toOutput ==1 && defined $baseFrac ){
        $toOutput = 0 if &fqParser::getBaseFraction($seq, $base) > $baseFrac;         
    }
    if($toOutput ==1 && defined $qualCutoff && defined $qualFrac ){
        $toOutput = 0 if &fqParser::getLowQualBaseFraction($qualSeq, $platform, $qualCutoff) > $qualFrac;
    }
    if($toOutput ==1 && defined $minAveQual){
        $toOutput = 0 if &fqParser::countAveQual($qualSeq, $platform) < $minAveQual;
    }
    if($toOutput ==1 && defined $tail){
        if( &fqParser::isBadTail($tail, $qualSeq, $platform)){
            $badTailN++;
            $toOutput = 0;
        }
    }
    if($toOutput == 1){
        print $readName;
        say $seq;
        print $plus;
        say $qualSeq;
        $passedN++;
    }elsif(defined $discardFile){
        print DIS $readName;
        say DIS $seq;
        print DIS $plus;
        say DIS $qualSeq;
    }
}
say STDERR "Total reads = $totalN";
say STDERR "Passed reads = "."$passedN (". sprintf ("%.2f", ($passedN/$totalN*100)) . "%)";
say STDERR "Bad tail reads = $badTailN" if defined $tail;

