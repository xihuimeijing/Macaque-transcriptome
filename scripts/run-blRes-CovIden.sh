#!/bin/sh
#######  bash run.sh refseqHuman+ 75 &
name=$1
identity=$2
cat ../blastp/result/refseq/blastp.human.out ../blastp/result/refseq/blastp.rhesus.out ../blastp/result/genbank/blastp.genbank.out > ${name}.out
blastpRes=$name.out
genomeTranGpe=/mnt/share/lisx/pacbio/NCBI_inputdata/protein-ORF/blastp/transcriptome_all/genome.transcriptome.gpe

perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f1 <(less -S ${genomeTranGpe} |cut -f 1,11,16,17|awk '{print $1"\t"$2"\t"($4-$3)/3-1}') -f2 <(less -S ${blastpRes} |grep -v "^#"|awk '{print $2"\t"$0}'|cut -f 1,2,4- ) -n 1 -e1 3 -e2 15|awk '{if ($6>20)print $0}'|cut -f 1-3,5-|awk '{print ($12-$11+1)/$3"\t"$0}' > ${name}.res
 
perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f <(less -S ${name}.res |awk '{if ($6>"'$identity'") print $0}'|cut -f 1,2,3,5|awk '{print $3"***"$2"***"$4"\t"$1}' )|sed 's/,/\t/g'|awk '{max=$2; if (NF>2){for (i=2;i<=NF;i++){ if ($i>=max ){max=$i}} print $1"\t"max} }' > tmp
perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f <(less -S ${name}.res |awk '{if ($6>"'$identity'") print $0}'|cut -f 1,2,3,5|awk '{print $3"***"$2"***"$4"\t"$1}' )|sed 's/,/\t/g'|awk '{if (NF==2){print $0} }'>>tmp
less -S tmp |sed 's/\*\*\*/\t/g' > ${name}.iden${identity}.allres
rm tmp
perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f <(less -S ${name}.iden${identity}.allres |awk '{print $1"***"$2"\t"$4}' ) |sed 's/,/\t/g'|awk '{max=$2; if (NF>2){for (i=2;i<=NF;i++){ if ($i>=max ){max=$i}} print $1"***"max} }' > ${name}.iden${identity}.topres
perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f <(perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f1 ${name}.iden${identity}.topres -f2 <(less -S ${name}.res|awk '{print $3"***"$2"***"$1"\t"$5","$12"-"$13"-"$4}' ) -n 1 -e1 0 -e2 1 |cut -f 1,3)|sed 's/,/\t/g'|awk '{printf $1"\t";for (i=1;i<=(NF-1)/2;i++){p=i*2;printf $p","};printf "\t";for (i=1;i<=(NF-1)/2;i++){p=i*2+1;printf $p","}; printf "\n" }'| sed 's/\*\*\*/\t/g'> ${name}.iden${identity}.Topres
#rm ${name}.iden${identity}.topres
#rm ${name}.iden${identity}.allres
#rm $name.out

#############

#cat refseq*Topres genbank.*Topres > all.Topres
#perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f <( perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f1 <(perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f <(less -S all.Topres |awk '{print $1"***"$2"\t"$3 }')|sed 's/,/\t/g'|awk '{max=$2; if (NF>2){for (i=2;i<=NF;i++){ if ($i>=max ){max=$i}}; print $1"***"max}; if (NF==2){ print $1"***"$2 } }'|sed 's/\*\*\*/\t/g' ) -f2 all.Topres -n 3 -e1 0 -e2 2 |cut -f 4-|awk '{print $1"***"$2"***"$3"\t"$4}')|sed 's/,,/,/g'|sed 's/\*\*\*/\t/g' > tmp 
#mv tmp all.Topres


#bash run-statistic.sh ${name}.iden${identity}.Topres
 
#perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f2 <(awk '{print $2"\t"$0}' ${name}.iden${identity}.Topres ) -f1 <(less -S /rd1/user/liym/transcriptome/visualization/geneNameAssign/version2_all/geneNameAssign.summary.tsv |awk '{print $1"\t"$7}'|sed 's/_/#/g'|sed 's/\./*/g'  ) -n 1 -e1 1 -e2 10|cut -f 2,5-|awk '{if ($2~/Single/ || $2~/L1/)print $0}' > L1Single.iden${identity}.Topres
#bash run-statistic.sh L1Single.iden${identity}.Topres

