#!/bin/sh
tissue=$1
gpe=${tissue}.cds.gpe
cage=../SR.cage.bed4
Junctionfile=../Junctionreads.RNAseq/newMonkey.all.junction.bed12
grep + $gpe > ${gpe}+
grep + $cage > ${cage}+
paste <(less -S ${gpe}+ |cut -f 1-5) <(cut -f 9 ${gpe}+ |sed 's/,/\t/g'|cut -f 1 ) <(cut -f 10 ${gpe}+ |sed 's/,/\t/g'|cut -f 1 ) <(cut -f 12 ${gpe}+ ) |awk '{print $2"\t"$6"\t"$7"\t"$3"\ttrans:"$3":"$4"-"$5"\t"$8"\t"$1}' > ${gpe}.all+exon.trans.gene
less -S ${gpe}.all+exon.trans.gene |awk '{if ($2<1000){$2=0}else{$2=$2-1000};print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > ${gpe}.firstexon+.1000
intersectBed -a ${cage}+ -b ${gpe}.firstexon+.1000 -wa -wb |awk '{if ($4==$8)print $0}'|awk '{print $1"\t"$2"\t"$6+1000"\tcage:"$4":"$2"-"$3"\t"$0}'|cut -f 1-4,13- > tmp.${gpe}.tmp
paste <(cut -f 1,2 tmp.${gpe}.tmp) <(less -S tmp.${gpe}.tmp|cut -f 5|sed 's/-/\t/g'|sed 's/:/\t/g'|cut -f 3 ) <(cut -f 4- tmp.${gpe}.tmp) >${gpe}.RawDistance.cage.trans+.gene
rm tmp.${gpe}.tmp
perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f  <( less -S ${gpe}.RawDistance.cage.trans+.gene |awk '{if ($2>$3 || $2==$3){print $6"\t"$7"\t"$4} }' |awk '{print $2"###"$1"\t"$3}' ) |awk '{print $1"\tlevel1\t"$2}'| sed 's/###/\t/g'|sort|uniq > ${gpe}.trans.gene.level1.cage.+
perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f1 ${gpe}.trans.gene.level1.cage.+ -f2 <(less -S ${gpe}.RawDistance.cage.trans+.gene |awk '{print $7"\t"$0}' ) -n 1 -e1 0 -e2 9 -only onlyf2|cut -f 2- > ${gpe}.RawDistance.level2.cage.trans+.gene

less -S $Junctionfile |cut -f 1,2,3,11,12|sed 's/,/\t/g'|awk '{print $1"\t"$2"\t"$2+$4"\t"$2+$7"\t"$3"\t"NR}'|cut -f 1,2,3,6|awk '{print $0"\tExon1"}' > ${gpe}.junctionreads
less -S $Junctionfile |cut -f 1,2,3,11,12|sed 's/,/\t/g'|awk '{print $1"\t"$2"\t"$2+$4"\t"$2+$7"\t"$3"\t"NR}'|cut -f 1,3,4,6|awk '{print $0"\tIntron"}' >> ${gpe}.junctionreads
less -S $Junctionfile |cut -f 1,2,3,11,12|sed 's/,/\t/g'|awk '{print $1"\t"$2"\t"$2+$4"\t"$2+$7"\t"$3"\t"NR}'|cut -f 1,4,5,6|awk '{print $0"\tExon2"}' >> ${gpe}.junctionreads
perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f <(perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f <(intersectBed -a ${gpe}.junctionreads -b ${gpe}.RawDistance.level2.cage.trans+.gene -wa -wb|cut -f 4,5,9,11,12|awk '{print $5"###"$4","$1","$3"\t"$2}'|sort|uniq) |awk '{if ($2=="Exon1,"|| $2=="Exon1,Exon2,Intron," || $2=="Exon1,Intron," || $2=="Intron," ) print $0}'|cut -f 1 |sed 's/,/\t/g'|cut -f 1,3|sort|uniq ) >${gpe}.level2.noJunction.transGene.cage+
 rm ${gpe}.junctionreads

less -S ${gpe}.level2.noJunction.transGene.cage+ |sed 's/,/\t/g'|sed 's/cage:+://g'|sed 's/-/\t/g'|awk '{a=2; for (i=3;i<=NF;i++) {if($i>$a)a=i };print $1"\tcage:+:"$a-1"-"$a  }'|sed 's/###/\t/g' > ${gpe}.nojunction.trans.gene.level2.cage.+ 

perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f2 <(less -S ${gpe}.RawDistance.level2.cage.trans+.gene |awk '{print $7"\t"$6"\t"$4"\t"$0}') -f1 ${gpe}.nojunction.trans.gene.level2.cage.+ -n 3 -e2 7 -e1 0|cut -f 7-|sort|uniq > tmp.level2.+
cut -f 1-3 tmp.level2.+ > distance+.bed3
bwtool summary distance+.bed3 /rd1/user/lisx/annotation/RNAseq/newMonkey/R-testis_HL7VFCCXX_L1/uniq.sorted.bw distance+.bed3.tmp
 
paste <(less -S tmp.level2.+ |awk '{print $7"\t"$6"\tlevel2\t"$4}') <(less -S distance+.bed3.tmp |awk '{print $1"\t"$2"\t"$3"\t"$5/$4}'|cut -f 4)|awk '{if ($5>0.7)print $1"\t"$2"\t"$3"\t"$4","$5}'|sort|uniq > ${gpe}.trans.gene.level2.cage.+



rm ${gpe}.firstexon+.1000
rm ${gpe}.RawDistance.*
rm ${gpe}.all+exon.trans.gene
rm ${gpe}.level2.noJunction.transGene.cage+
rm ${gpe}+ 
rm ${gpe}.nojunction.trans.gene.level2.cage.+
rm distance+.*
rm tmp.level2.+






