#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
my ($infilea,$infileb,$infilec,$n,$splice);
my (@line,%ref,%pacbioread,%inter,%output);
my ($ref,$pacbio);

GetOptions(
	   'f1|infilea=s'       => \$infilea,
	  ) || usage();


 #### the refgene gpe file
open INDEX,"$infilea"or die "Can't open file $infilea:$!";
while (<INDEX>){
    chomp;
    @line=split (/\t/,$_);
    $ref=$_;
    my $reflength=gpeLength($_); 
    my $length5UTR=gpe5UTRlenth($_);
    my $length3UTR=gpe3UTRlenth($_);
    if ($line[3] eq "-"){ my $tmp=$length3UTR; $length3UTR=$length5UTR;$length5UTR=$tmp }
    if ($line[6]==$line[7]){
	say "0:0:0:",$reflength,"\t",$_;
    }elsif($line[5]==$line[7] && $line[4]==$line[6]){
	say "0:",$reflength,":0:",$reflength,"\t",$_;
    }else{
	say $length5UTR,":",$reflength-$length5UTR-$length3UTR,":",$length3UTR,":",$reflength,"\t",$_;
    }
}
sub gpeLength{
    my @ref=split(/\t/,$_[0]); 
    my $n_ref=$ref[8];
    my $transStart; my $transEnd; my $CDSstart; my $CDSend; my @exonStart; my @exonEnd;
    $transStart=$ref[4]; $transEnd=$ref[5]; $CDSstart=$ref[6]; $CDSend=$ref[7]; @exonStart=split (/,/,$ref[9]); @exonEnd=split(/,/,$ref[10]);
    my $length=0;
    for (my $i=0;$i<$n_ref;$i++){
	my $exon=$exonEnd[$i]-$exonStart[$i];
	$length=$length+$exon;
    }
    return $length;
}
sub gpe5UTRlenth{
    my @ref=split(/\t/,$_[0]); 
    my $n_ref=$ref[8];
    my $cds;
    my $transStart; my $transEnd; my $CDSstart; my $CDSend; my @exonStart; my @exonEnd;
    $transStart=$ref[4]; $transEnd=$ref[5]; $CDSstart=$ref[6]; $CDSend=$ref[7]; @exonStart=split (/,/,$ref[9]); @exonEnd=split(/,/,$ref[10]);
    for (my $i=0;$i<$n_ref;$i++){
	if ($exonStart[$i]<=$CDSstart && $exonEnd[$i]>=$CDSstart){ $cds=$i}
    }
   my @UTRstart;my @UTRend;
    for (my $i=0;$i<$cds+1;$i++){
	$UTRend[$i]=$exonEnd[$i];
	$UTRstart[$i]=$exonStart[$i];
	if ($i==$cds){$UTRend[$i]=$CDSstart}
    }
    my $length=0;
    for (my $i=0;$i<$cds+1;$i++){ $length=$length+$UTRend[$i]-$UTRstart[$i]}
    return $length;
}

sub gpe3UTRlenth{
    my @ref=split(/\t/,$_[0]); 
    my $n_ref=$ref[8];
    my $cds;
    my $transStart; my $transEnd; my $CDSstart; my $CDSend; my @exonStart; my @exonEnd;
    $transStart=$ref[4]; $transEnd=$ref[5]; $CDSstart=$ref[6]; $CDSend=$ref[7]; @exonStart=split (/,/,$ref[9]); @exonEnd=split(/,/,$ref[10]);
    for (my $i=0;$i<$n_ref;$i++){
	if ($exonStart[$i]<=$CDSend && $exonEnd[$i]>=$CDSend){ $cds=$i}
    }
   my @UTRstart;my @UTRend;
    for (my $i=$cds;$i<$n_ref;$i++){
	$UTRend[$i]=$exonEnd[$i];
	$UTRstart[$i]=$exonStart[$i];
	if ($i==$cds){$UTRstart[$i]=$CDSend}
    }
    my $length=0;
    for (my $i=$cds;$i<$n_ref;$i++){ $length=$length+$UTRend[$i]-$UTRstart[$i]}
    if ($CDSend==$CDSstart){$length=0; for (my $i=0;$i<$n_ref;$i++){ $length=$length+$exonEnd[$i]-$exonStart[$i]}	}
    return $length;
}
sub gpe5UTRstart{
    my @ref=split(/\t/,$_[0]); 
    my $n_ref=$ref[8];
    my $cds;
    my $transStart; my $transEnd; my $CDSstart; my $CDSend; my @exonStart; my @exonEnd;
    $transStart=$ref[4]; $transEnd=$ref[5]; $CDSstart=$ref[6]; $CDSend=$ref[7]; @exonStart=split (/,/,$ref[9]); @exonEnd=split(/,/,$ref[10]);
    for (my $i=0;$i<$n_ref;$i++){
	if ($exonStart[$i]<=$CDSstart && $exonEnd[$i]>=$CDSstart){ $cds=$i}
    }
   my @UTRstart;my @UTRend;
    for (my $i=0;$i<$cds+1;$i++){
	$UTRend[$i]=$exonEnd[$i];
	$UTRstart[$i]=$exonStart[$i];
	if ($i==$cds){$UTRend[$i]=$CDSstart}
    }
    return @UTRstart;
}
sub gpe5UTRend{
    my @ref=split(/\t/,$_[0]); 
    my $n_ref=$ref[8];
    my $cds;
    my $transStart; my $transEnd; my $CDSstart; my $CDSend; my @exonStart; my @exonEnd;
    $transStart=$ref[4]; $transEnd=$ref[5]; $CDSstart=$ref[6]; $CDSend=$ref[7]; @exonStart=split (/,/,$ref[9]); @exonEnd=split(/,/,$ref[10]);
    for (my $i=0;$i<$n_ref;$i++){
	if ($exonStart[$i]<=$CDSstart && $exonEnd[$i]>=$CDSstart){ $cds=$i}
    }
   my @UTRstart;my @UTRend;
    for (my $i=0;$i<$cds+1;$i++){
	$UTRend[$i]=$exonEnd[$i];
	$UTRstart[$i]=$exonStart[$i];
	if ($i==$cds){$UTRend[$i]=$CDSstart}
    }
    return @UTRend;
}
sub gpe3UTRstart{
    my @ref=split(/\t/,$_[0]); 
    my $n_ref=$ref[8];
    my $cds;
    my $transStart; my $transEnd; my $CDSstart; my $CDSend; my @exonStart; my @exonEnd;
    $transStart=$ref[4]; $transEnd=$ref[5]; $CDSstart=$ref[6]; $CDSend=$ref[7]; @exonStart=split (/,/,$ref[9]); @exonEnd=split(/,/,$ref[10]);
    for (my $i=0;$i<$n_ref;$i++){
	if ($exonStart[$i]<=$CDSend && $exonEnd[$i]>=$CDSend){ $cds=$i}
    }
   my @UTRstart;my @UTRend; my $l=0;my @outstart;my @outend;
    for (my $i=$cds;$i<$n_ref;$i++){
	$UTRend[$i]=$exonEnd[$i];
	$UTRstart[$i]=$exonStart[$i];
	if ($i==$cds){$UTRstart[$i]=$CDSend}
	$outstart[$l]=$UTRstart[$i];$outend[$l]=$UTRend[$i];$l++;
    }
    return @outstart;
}
sub gpe3UTRend{
    my @ref=split(/\t/,$_[0]); 
    my $n_ref=$ref[8];
    my $cds;
    my $transStart; my $transEnd; my $CDSstart; my $CDSend; my @exonStart; my @exonEnd;
    $transStart=$ref[4]; $transEnd=$ref[5]; $CDSstart=$ref[6]; $CDSend=$ref[7]; @exonStart=split (/,/,$ref[9]); @exonEnd=split(/,/,$ref[10]);
    for (my $i=0;$i<$n_ref;$i++){
	if ($exonStart[$i]<=$CDSend && $exonEnd[$i]>=$CDSend){ $cds=$i}
    }
   my @UTRstart;my @UTRend;my $l=0;my @outstart;my @outend;
    for (my $i=$cds;$i<$n_ref;$i++){
	$UTRend[$i]=$exonEnd[$i];
	$UTRstart[$i]=$exonStart[$i];
	if ($i==$cds){$UTRstart[$i]=$CDSend}
	$outstart[$l]=$UTRstart[$i];$outend[$l]=$UTRend[$i];$l++;
    }
    return @outend;
}
sub gpeStart{
    my @ref=split(/\t/,$_[0]); 
    my $n_ref=$ref[8];
    my $transStart; my $transEnd; my $CDSstart; my $CDSend; my @exonStart; my @exonEnd;
    $transStart=$ref[4]; $transEnd=$ref[5]; $CDSstart=$ref[6]; $CDSend=$ref[7]; @exonStart=split (/,/,$ref[9]); @exonEnd=split(/,/,$ref[10]);
    return @exonStart;
}
sub gpeEnd{
    my @ref=split(/\t/,$_[0]); 
    my $n_ref=$ref[8];
    my $transStart; my $transEnd; my $CDSstart; my $CDSend; my @exonStart; my @exonEnd;
    $transStart=$ref[4]; $transEnd=$ref[5]; $CDSstart=$ref[6]; $CDSend=$ref[7]; @exonStart=split (/,/,$ref[9]); @exonEnd=split(/,/,$ref[10]);
    return @exonEnd;
}
sub gpeCDSstart{
    my @ref=split(/\t/,$_[0]); 
    my $n_ref=$ref[8];
    my $transStart; my $transEnd; my $CDSstart; my $CDSend; my @exonStart; my @exonEnd; my @CDSst; my @CDSen; 
    $transStart=$ref[4]; $transEnd=$ref[5]; $CDSstart=$ref[6]; $CDSend=$ref[7]; $n=$ref[8];@exonStart=split (/,/,$ref[9]); @exonEnd=split(/,/,$ref[10]);
    my ($a,$b);
    for (my $i=0; $i<$n;$i++){
	if ($exonStart[$i]<=$CDSstart && $exonEnd[$i]>=$CDSstart){ $a=$i}
	if ($exonStart[$i]<=$CDSend && $exonEnd[$i]>=$CDSend){ $b=$i}
    }
    my $l=0;
    for (my $i=$a;$i<=$b;$i++){
	if ($i==$a){ $CDSst[$l]=$CDSstart }else{$CDSst[$l]=$exonStart[$i]}
	if ($i==$b){ $CDSen[$l]=$CDSend }else{ $CDSen[$l]=$exonEnd[$i] }
	$l++
    }
    return @CDSst;
}
sub gpeCDSend{
    my @ref=split(/\t/,$_[0]); 
    my $n_ref=$ref[8];
    my $transStart; my $transEnd; my $CDSstart; my $CDSend; my @exonStart; my @exonEnd; my @CDSst; my @CDSen; 
    $transStart=$ref[4]; $transEnd=$ref[5]; $CDSstart=$ref[6]; $CDSend=$ref[7]; $n=$ref[8];@exonStart=split (/,/,$ref[9]); @exonEnd=split(/,/,$ref[10]);
    my ($a,$b);
    for (my $i=0; $i<$n;$i++){
	if ($exonStart[$i]<=$CDSstart && $exonEnd[$i]>=$CDSstart){ $a=$i}
	if ($exonStart[$i]<=$CDSend && $exonEnd[$i]>=$CDSend){ $b=$i}
    }
    my $l=0;
    for (my $i=$a;$i<=$b;$i++){
	if ($i==$a){ $CDSst[$l]=$CDSstart }else{$CDSst[$l]=$exonStart[$i]}
	if ($i==$b){ $CDSen[$l]=$CDSend }else{ $CDSen[$l]=$exonEnd[$i] }
	$l++
    }
    return @CDSen;
}

sub bedStart{
    my @ref=split(/\t/,$_[0]); 
    my $n_ref=$ref[9];
    my $transStart; my $transEnd; my $CDSstart; my $CDSend; my @exonStart; my @exonEnd;
    $transStart=$ref[1]; $transEnd=$ref[2]; $CDSstart=$ref[6]; $CDSend=$ref[7]; @exonStart=split (/,/,$ref[9]); @exonEnd=split(/,/,$ref[10]);
	my @start=split(/,/,$ref[11]); my @end=split(/,/,$ref[10]);
    for (my $i=0;$i<$n_ref;$i++){ $exonStart[$i]=$transStart+$start[$i]; $exonEnd[$i]=$exonStart[$i]+$exonEnd[$i] }
    return @exonStart;
}
sub bedEnd{
    my @ref=split(/\t/,$_[0]); 
    my $n_ref=$ref[9];
    my $transStart; my $transEnd; my $CDSstart; my $CDSend; my @exonStart; my @exonEnd;
    $transStart=$ref[1]; $transEnd=$ref[2]; $CDSstart=$ref[6]; $CDSend=$ref[7]; @exonStart=split (/,/,$ref[9]); @exonEnd=split(/,/,$ref[10]);
	my @start=split(/,/,$ref[11]); my @end=split(/,/,$ref[10]);
    for (my $i=0;$i<$n_ref;$i++){ $exonStart[$i]=$transStart+$start[$i]; $exonEnd[$i]=$exonStart[$i]+$exonEnd[$i] }
    return @exonEnd;
}
sub usage{
print <<HELP;
Usage:  perl $0 -f file >file
Author:	Shuxian Li 2017-09-28
Options:
	-f1	File	the refgene.gpe(with bin)
	-h    --help	Print this help information.
Output: 5UTR:CDS:3UTR:length gpe 
HELP
    exit(-1);
}