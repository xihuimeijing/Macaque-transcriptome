#!/bin/sh
xlab=testis
ylab=brain
folder=compare
cat head <(perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f1 ${ylab}/gene.3utr -f2 ${xlab}/gene.3utr -n 1 -e1 1 -e2 1 |cut -f 1,2,4|awk '{if (($2-$3)>-30 && ($2-$3)<30 )print $0}'|cut -f 2-|awk '{print log($1)/log(2)"\t"log($2)/log(2)}') > ${folder}/${ylab}-${xlab}.grey
cat head <(perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f1 ${ylab}/gene.3utr -f2 ${xlab}/gene.3utr -n 1 -e1 1 -e2 1 |cut -f 1,2,4|awk '{if (($2-$3)<-30 )print $0}' |cut -f 2-|awk '{print log($1)/log(2)"\t"log($2)/log(2)}' ) > ${folder}/${ylab}-${xlab}.blue
cat head <(perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f1 ${ylab}/gene.3utr -f2 ${xlab}/gene.3utr -n 1 -e1 1 -e2 1 |cut -f 1,2,4|awk '{if (($2-$3)>30 )print $0}'|cut -f 2-|awk '{print log($1)/log(2)"\t"log($2)/log(2)}' ) > ${folder}/${ylab}-${xlab}.red
 
#Rscript Rplot.R ${folder}/${ylab}-${xlab}.grey ${folder}/${ylab}-${xlab}.blue ${folder}/${ylab}-${xlab}.red $xlab $ylab ${ylab}-${xlab}.pdf





