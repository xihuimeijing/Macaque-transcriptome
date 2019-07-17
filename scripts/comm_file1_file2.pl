#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
my ($infilea,$infileb,$n,$extend1,$extend2,$only);
my (@line,%filea,%fileb);

GetOptions(
	   'f1|infilea=s'       => \$infilea,
	   'f2|infileb=s'       => \$infileb,
	   'n=i'               => \$n,
	   'e1=i'               => \$extend1,
	   'e2=i'		=> \$extend2,
	   'only=s'		=> \$only,
           'h|help'          => sub{usage()}
	  ) || usage();

if ((defined $only)==0){ $only="comm"}
if ($only eq "onlyf1"){
    my $tmp=$infilea;$infilea=$infileb;$infileb=$tmp; 
    $tmp=$extend1;$extend1=$extend2;$extend2=$tmp;    
}
open INDEX,"$infilea"or die "Can't open file $infilea:$!";
while (<INDEX>){
    chomp;
    @line=split (/\t/,$_);
    my $nrow;
    my $ex;
    if ($n==1){
	$nrow=$line[0];
    }else{
	for (my $i=0;$i<$n;$i++){
	    if ($i==$n-1){
		$nrow.=$line[$i];
	    }elsif ($i==0){
		$nrow=$line[$i]."\t";
	    }else{
		$nrow=$nrow.$line[$i]."\t";
	    }
	}
    }
    
    if ($extend1!=0){
	my $a=@line;
	for (my $i=$n;$i<$a;$i++){
	    if ($i==$a-1){
		$ex.=$line[$i];
	    }elsif ($i==$n){
		$ex=$line[$i]."\t";
	    }else{
		$ex=$ex.$line[$i]."\t";
	    }
	}
	$filea{$nrow}=$ex;
    }elsif($extend1==0){
	$filea{$nrow}=undef;
    }
}

open INDEX,"$infileb"or die "Can't open file $infileb:$!";
while (<INDEX>){
    chomp;
    @line=split (/\t/,$_);
    my $nrow=undef;
    if ($n==1){
	$nrow=$line[0];
    }else{
	for (my $i=0;$i<$n;$i++){
	    if ($i==$n-1){
		$nrow.=$line[$i];
	    }elsif ($i==0){
		$nrow=$line[$i]."\t";
	    }else{
		$nrow=$nrow.$line[$i]."\t";
	    }
	}
    }
    if (exists $filea{$nrow} ){
	if ($only eq "comm"){
	    if ($extend1!=0 ){
		say $nrow,"\t",$filea{$nrow},"\t",$_;
	    }else{
	    	say $nrow,"\t",$_;
	    }
	}
    }else{
	if ($only ne "comm"){ say $_}
    }
}

sub usage{
print <<HELP;
Usage:  perl $0 -f1 file1 -f2 file2 -n number -e1 number -e1 number >file
Author:	Shuxian Li 2015-11-13
Options:
	-f1	File	the file that contain n row
	-f2	File	the file that contain n+ row
	-n	Int	the same rows which the two file contains
	-e1	Int	the extend number of rows which the first file (file1) contains
	-e2	Int	the extend number of rows which the first file (file2) contains
	-only   string 	comm or onlyf1 or onlyf2 : output lines is comm or only in file1 or file2
	-h    --help	Print this help information.
Output: the comm file  
HELP
    exit(-1);
}
