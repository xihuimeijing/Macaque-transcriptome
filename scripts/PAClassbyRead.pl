#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my $assignRead;
GetOptions(
	    'a|assign=s'          => \$assignRead,
            'h|help'              => sub{usage()}
          ) || usage();
open IN, "$ARGV[0]" or die "Can't open file $ARGV[0]: $!";
open BED, "$assignRead" or die "Can't open file $assignRead: $!";


sub usage{
print <<HELP;
Usage: perl scriptName -a assignRead PA.bed > PAClassify.bed6+
Output:chr, PAstart, PAend, MEreadID,M MEreadNO,strand ...,type, transID
Options:
   -a --assign         Assign reads in bed6+
   -h --help           Show this help information
HELP
   exit(-1);
}

my(%ReadConExon,%ReadTran,%ReadExon,%Read);
while(<BED>){
   chomp;
   my @filed=split;
   my ($chr,$readStart,$readEnd,$readID,$strand,$exonNum,$blockSize,$blockStart,$transID,$conExonNum,$coverage,$identify)=@filed[0..3,5,9,10,11,15,17,18,19];
   $ReadConExon{$readID}=$conExonNum;
   $ReadTran{$readID}=$transID;
   $ReadExon{$readID}=$exonNum;
   $Read{$readID}=$_;
}

my (%PA,%PAType,%PAME,%PASE,%PATrans,$key,$value,$type);
while(<IN>){
  chomp;
  my @filed=split;
  my($chr,$Start,$End,$read,$count,$strand)=@filed[0..5];
  my @reads=split(/,/,$read);
  $key=join("\t",$chr,$Start,$End,$strand);
  $value=join("\t",@filed[3..$#filed]);

  for(my $i=0;$i<@reads;$i++){
     if($ReadExon{$reads[$i]}){
        if($ReadExon{$reads[$i]}>1){
          $type="ME";
          $PAME{$key}.="$reads[$i],";
        }else{
          $type="SE";
          $PASE{$key}.="$reads[$i],";
        }
        $PATrans{$key}.="$ReadTran{$reads[$i]},";
     }
     $PAType{$key}.="$type,";
   }
   $PA{$key}=$value;
}
my (@key,@Utrans,%trans,@PAME,@PASE,@tmp);
foreach my $PA(keys %PAType){
   @key=split(/\t/,$PA);
   @tmp=split(/\t/,$PA{$PA});
   my @trans=split(/,/,$PATrans{$PA});
   undef %trans;
   undef @Utrans;
   undef @PAME;
   undef @PASE;
   for(my $i=0;$i<@trans;$i++){
      if($trans{$trans[$i]}){
         next;
      }else{
        $trans{$trans[$i]}=1;
        push @Utrans,$trans[$i];
      }
   }

      if($PAType{$PA}=~/SE/ && $PAType{$PA}!~/ME/){
         @PASE=split(/,/,$PASE{$PA});
         print join("\t",@key[0..2],join(",",@PASE),$#PASE+1,@tmp[2..$#tmp],"SE",join(",",@Utrans),"\n");
      }
      if($PAType{$PA}=~/SE/ && $PAType{$PA}=~/ME/){
        @PAME=split(/,/,$PAME{$PA});
        print join("\t",@key[0..2],join(",",@PAME),$#PAME+1,@tmp[2..$#tmp],"ME&SE",join(",",@Utrans),"\n");
        @PASE=split(/,/,$PASE{$PA});
        for(my $j=0; $j<@PASE; $j++){
          if($Read{$PASE[$j]}){
            print STDERR join ("\t",$Read{$PASE[$j]},$key[0].":".$key[1]."-".$key[2],$#PASE+1,"\n");
          }
        }
      }
      if($PAType{$PA}=~/ME/ && $PAType{$PA}!~/SE/){
         @PAME=split(/,/,$PAME{$PA});
          print join("\t",@key[0..2],join(",",@PAME),$#PAME+1,@tmp[2..$#tmp],"ME",join(",",@Utrans),"\n");
     }
}

