#!/usr/bin/env perl
use 5.010;
use warnings;
use Getopt::Long;
use File::Basename;
use List::MoreUtils qw(uniq);

my ($gpeFile, $bin);
GetOptions(
            'g|gpe=s'   => \$gpeFile,
            'b|bin'     => \$bin,
            'h|help'    => sub{usage()}
        ) || usage();
$ARGV[0] = '-' unless defined $ARGV[0];
open IN, "$ARGV[0]" or die "Can't read file ($ARGV[0]): $!";
open GPE, $gpeFile or die "Can't read file ($gpeFile): $!";
my %transExonNum;
while(<GPE>){
	chomp;
	my @split=split;
	if(defined $bin){
		$transExonNum{$split[1]}=$split[8];
	}else{
		$transExonNum{$split[0]}=$split[7];
	}
}
my $novelG=0;
while(<IN>){
	chomp;
	my @split=split;
	my($blockNum,$tag,$transs,$genes,$cssIntronN,$covOnRead,$covOnTrans)=@split[9,14..19];
	my @transNames=split /,/,$transs;
	my @geneNames=split /,/,$genes;
	my @cssIntronN=split /,/,$cssIntronN;
	my @covOnTrans=split /,/,$covOnTrans;
	my @covOnRead=split /,/,$covOnRead;
	my @uniqGeneNames=uniq @geneNames;
	if($#uniqGeneNames>0){
		say STDERR join "\t",$_,(join ",",@uniqGeneNames);
	}else{
		if($tag eq "IG"){
			$novelG+=1;
			my $geneName="Novel-".$novelG;
			say "$_\t$geneName";
		}elsif($tag eq "E"){
			if($blockNum==1){ #SE reads
				my $unambi=1;
				for(my $i=0;$i<=$#transNames;$i++){
					if($transExonNum{$transNames[$i]}==1 && $covOnRead[$i]>=0.5){
						$unambi=0;
						say join "\t",$_,$uniqGeneNames[0];
						last;
					}
				}
				if($unambi==1){
					say STDERR "$_\t$uniqGeneNames[0]";
				}
			}else{ #ME reads
				say "$_\t$uniqGeneNames[0]";
			}
		}elsif($tag eq "I"){
			if($blockNum==1){ #SE reads
				my $unambi=1;
				for(my $i=0;$i<=$#transNames;$i++){
					if($transExonNum{$transNames[$i]}==1 && $covOnRead[$i]>=0.5){
						$unambi=0;
						say join "\t",$_,$uniqGeneNames[0];
						last;
					}
				}
				if($unambi==1){
					say STDERR "$_\t$uniqGeneNames[0]";
				}
			}else{ #ME reads
				my $unambi=1;
				for(my $i=0;$i<=$#transNames;$i++){
					if($transExonNum{$transNames[$i]}>1 && $covOnRead[$i]>=0.8){
						$unambi=0;
						say join "\t",$_,$uniqGeneNames[0];
						last;
					}
				}
				if($unambi==1){
					say STDERR "$_\t$uniqGeneNames[0]";
				}
			}
		}
	}
}
sub usage{
    my $scriptName = basename $0;
print <<HELP;
Usage: perl $scriptName INPUT >OUTPUT
    INPUT is the output result of readsAssigner.pl. If INPUT isn't specified, input from STDIN.
Options:
    -g --gpe    FILE    The gene annotation file in gpe format
    -b --bin            With bin column in --gpe
    -h --help           Print this help information
HELP
    exit(-1);
}