#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use vars qw ($putative $gpe  $bin);
use vars qw ($fh1 $fh2);
GetOptions (
        "putativeCDS|i:s"      =>      \$putative,
        "gpe|g:s"              =>      \$gpe,
        "bin|b"                =>      \$bin,
        "help|h"               =>      sub{&usage;exit(0);}
);

if (! open $fh1,"< $putative"){die "cannot open $putative!";}
if (! open $fh2,"< $gpe"){die "cannot open $gpe!";}

my($lineA,@pos,@nameA,%hash,$lineB,@array,@cds,@exonstart,@exonend,$cdsstart,$cdsend,$m);
while ($lineA=<$fh1>)
    {
        @pos=split/\t/,$lineA;
        @nameA=split/\|/,$pos[0];
        $hash{$nameA[1]}=$pos[1]."|".$pos[2];   
    }

while ($lineB=<$fh2>)
    {
        
        @array=split/\t/,$lineB;
        if (defined $bin) {shift @array;}
        @cds=split/\|/,$hash{$array[0]};
        if ($cds[0] eq "-" && $cds[1] eq "-") {
        $cdsstart=$array[4];
        $cdsend=$array[4];
        }
        else
        {
        @exonstart=split/,/,$array[8];
        @exonend=split/,/,$array[9];
        if ($array[2] eq '+')
            {
                for (my $i=0;$i<$array[7];$i++)
                    {
                        if ($exonstart[$i]+$cds[0]-1<$exonend[$i])
                            {
                                $cdsstart=$exonstart[$i]+$cds[0]-1;
                                last;
                            }
                        else
                            {
                                $cds[0]=$cds[0]-$exonend[$i]+$exonstart[$i];    
                            }
                    }
                for (my $j=0;$j<$array[7];$j++)
                    {
                        if ($exonstart[$j]+$cds[1]-1<$exonend[$j])
                            {
                                $cdsend=$exonstart[$j]+$cds[1];
                                last;
                            }
                        else
                            {
                                $cds[1]=$cds[1]-$exonend[$j]+$exonstart[$j];    
                            }
                    }   
                    
            }
        elsif ($array[2] eq '-')
            {
                for (my $i=$array[7]-1;$i>=0;$i--)
                    {
                        if ($exonend[$i]-$cds[0]+1>$exonstart[$i])
                            {
                                $cdsend=$exonend[$i]-$cds[0]+1;
                                last;
                            }
                        else
                            {
                                $cds[0]=$cds[0]-$exonend[$i]+$exonstart[$i];       
                            }
                    }
                for (my $j=$array[7]-1;$j>=0;$j--)
                    {
                        if ($exonend[$j]-$cds[1]+1>$exonstart[$j])
                            {
                                $cdsstart=$exonend[$j]-$cds[1];
                                last;
                            }
                        else
                            {
                                $cds[1]=$cds[1]-$exonend[$j]+$exonstart[$j];       
                            }
                    }    
            }
        }
        $array[5]=$cdsstart;
        $array[6]=$cdsend;
        print join "\t",@array;
    }






sub usage{
print STDERR <<HELP
Usage: perl $0 -i * -g * -b 
        --putativeCDS|-i:s         putativeCDS file produced by putativeTranscript.pl; 
        --gpe|-g:s                 genePredExt;
        --bin|-b                   gpe with bin
        --help|-h                  print this message
HELP
} 
