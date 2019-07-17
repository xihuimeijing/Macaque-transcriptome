#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use Bio::Perl;
use Bio::DB::Fasta;
my($bed,$refGene,$upstream,$downstream,$strand);
my $opt=GetOptions(
                        'b|bed=s'       => \$bed,
                        'u|up=i'        => \$upstream,
                        'd|down=i'      => \$downstream,
                        'f|fa=s'        => \$refGene,
                        's|strand'      => \$strand,
                        'h|help'	=> sub{&usage;exit(-1);}
                  );
my ($BED,%hash,$process);
if(defined $bed){
    open $BED,"$bed" or die "Can't open file $bed:$!";
}else{
    $BED=\*STDIN;
}
my $db = Bio::DB::Fasta->new($refGene);
for(my $i=$upstream;$i<=$downstream;$i++){
    $hash{$i}{'A'}=0;
    $hash{$i}{'T'}=0;
    $hash{$i}{'G'}=0;
    $hash{$i}{'C'}=0;
}
while(<$BED>){
    chomp;
    my @split=split;
    my $chr=$split[0];
    my @ref_seq;
    if(defined $strand){
        if($split[5] eq '+'){
            my $start=$split[1]+$upstream+1; #1-base corrdinate
            my $end=$split[2]+$downstream;
            if($start>0){
                @ref_seq=split //,uc($db->seq($chr,$start,$end));
            }else{
                next;
            }
        }else{
            my $start=$split[1]-$downstream+1; #1-base corrdinate
            my $end=$split[2]-$upstream;
            if($start>0){
                @ref_seq=split //,uc($db->seq($chr,$end,$start));
            }else{
                next;
            }
        }
    }else{
        my $start=$split[1]+$upstream+1; #1-base corrdinate
        my $end=$split[2]+$downstream;
        if($start>0){
            @ref_seq=split //,uc($db->seq($chr,$start,$end));
        }else{
            next;
        }
    }
    my $tag=$upstream;
    for(my $i=0;$i<=$#ref_seq;$i++){
        if($ref_seq[$i] eq 'A'){
            $hash{$tag}{'A'}+=1;
        }elsif($ref_seq[$i] eq 'T'){
            $hash{$tag}{'T'}+=1;
        }elsif($ref_seq[$i] eq 'G'){
            $hash{$tag}{'G'}+=1;
        }elsif($ref_seq[$i] eq 'C'){
            $hash{$tag}{'C'}+=1;
        }
        $tag++;
    }
}
print "Position\t";
foreach my $key1(sort keys %{$hash{'0'}}){
    print "$key1\t";
}
print "\n";
foreach my $pos(sort keys %hash){
    print "$pos\t";
    foreach my $base(sort keys %{$hash{$pos}}){
        print "$hash{$pos}{$base}\t";
    }
    print "\n";
}

sub usage{
print STDERR <<HELP 
Usage:	perl $0 -b <*.bed> -f hg19.fa >result.file 2>log
        Statistics the nucleotides count for the input point and upstream and downstream regions. 
Output: relative_position base_count(4 columns, A, C, G, T) 
        'b|bed'    FILE     The input file in bed3 or bed3+ format.(If not given, it can be read from STDIN)
        'u|up'     INT      The upstream region, negative value
        'd|down'   INT      The downstream region, positive value
        'f|fa'     FILE     Reference genome, fasta format
        's|strand' LOGIC    Considering strandness information,this restricts the input data in bed6/bed6+ format
        'help|h'            Print this help message    
HELP
}