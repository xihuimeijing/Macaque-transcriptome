#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min);
use vars qw ($input $output);
use vars qw ($fh1 $fh2 );
GetOptions (
	"Fasta|i:s"             =>      \$input,
	"output|o:s"            =>      \$output,
	"help|h"                =>      sub{&usage;exit(0);}
);

open $fh1,"< $input";
open $fh2,"> $output";

my($lineA,@seq,$codon,$aaNum,$BaseNum,$aa,$name);
while ($lineA=<$fh1>)
    {
    chomp $lineA;
    if ($lineA=~/>/){$name=$lineA;next;}
    else
        {
        @seq=split/|/,$lineA;
        $BaseNum=@seq;
        $aaNum=$BaseNum/3;
        for (my$i=0;$i<$aaNum;$i++)
            {
            $codon=$seq[3*$i];
            $codon.=$seq[3*$i+1];
            $codon.=$seq[3*$i+2];
            if ($codon=~/(GCT)|(GCC)|(GCA)|(GCG)/i) {$aa.='A';}
            elsif ($codon=~/(CGT)|(CGC)|(CGA)|(CGG)|(AGA)|(AGG)/i){$aa.='R';}
            elsif ($codon=~/(AAT)|(AAC)/i){$aa.='N';}
            elsif ($codon=~/(GAT)|(GAC)/i){$aa.='D';}
            elsif ($codon=~/(TGT)|(TGC)/i){$aa.='C';}
            elsif ($codon=~/(CAA)|(CAG)/i){$aa.='Q';}
            elsif ($codon=~/(GAA)|(GAG)/i){$aa.='E';}
            elsif ($codon=~/(GGT)|(GGC)|(GGA)|(GGG)/i){$aa.='G';}
            elsif ($codon=~/(CAT)|(CAC)/i){$aa.='H';}
            elsif ($codon=~/(ATT)|(ATC)|(ATA)/i){$aa.='I';}
            elsif ($codon=~/(TTA)|(TTG)|(CTT)|(CTC)|(CTA)|(CTG)/i){$aa.='L';}
            elsif ($codon=~/(AAA)|(AAG)/i){$aa.='K';}
            elsif ($codon=~/ATG/i){$aa.='M';}
            elsif ($codon=~/(TTT)|(TTC)/i){$aa.='F';}
            elsif ($codon=~/(CCT)|(CCC)|(CCA)|(CCG)/i){$aa.='P';}
            elsif ($codon=~/(TCT)|(TCC)|(TCA)|(TCG)|(AGT)|(AGC)/i){$aa.='S';}
            elsif ($codon=~/(ACT)|(ACC)|(ACA)|(ACG)/i){$aa.='T';}
            elsif ($codon=~/(TGG)/i){$aa.='W';}
            elsif ($codon=~/(TAT)|(TAC)/i){$aa.='Y';}
            elsif ($codon=~/(GTT)|(GTC)|(GTA)|(GTG)/i){$aa.='V';}
            elsif ($codon=~/(TAA)|(TGA)|(TAG)/i){$aa.='*';}
            }
        print $fh2 "$name\n";
        print $fh2 "$aa\n";
        $codon=undef;$name=undef;$aa=undef;@seq=undef;$BaseNum=undef;$aaNum=undef;
        }
        
    }
    
sub usage{
print STDERR <<HELP
Usage: perl $0 -i [.fa] -o [output] -h
        --Fasta|-i:s                     input;fasta format;     
        --output|-o:s                    output
        --help|-h                        print this message
HELP
}

#this script was used to convert base.fa format to protein.format, with 0-based start;
#2013.4.3
