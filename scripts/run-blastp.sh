#!/bin/sh
###
### extract aa sequence of bgm macaque trascriptome  
	### get base sequence 
bash run1-genome-makeblastdb.sh all
	### base sequence --> protein
perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f2 <(less -S transcriptome_all/all.transcriptome.fa |awk '{if(0==NR%2)printf("%s\n",$0);else printf("%s\t",$0)}' ) -f1 <(less -S transcriptome_all/genome.transcriptome.gpe |awk '{print ">"$13"\t>"$1}' ) -n 1 -e1 1 -e2 1 |cut -f 2,4|sed 's/>//g'> transcript.fa
perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f2 <(less -S transcript.fa ) -f1 <(less -S transcriptome_all/genome.transcriptome.gpe |awk '{if ($6!=$7)print $0}'|cut -f 1,3,15-17 ) -n 1 -e1 5 -e2 1 |cut -f 1-5,7|grep "+"|awk '{print $1"\t"substr($6,$4+1,$5-$4)}' > transcript.putativecds+.fa
perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f2 <(less -S transcript.fa ) -f1 <(less -S transcriptome_all/genome.transcriptome.gpe |awk '{if ($6!=$7)print $0}'|cut -f 1,3,15-17 ) -n 1 -e1 5 -e2 1 |cut -f 1-5,7|grep "-"|awk '{print $1"\t"substr($6,$4+1,$5-$4)}' > transcript.putativecds-.fa
bash run-test.sh |sed 'y/ATCG/atcg/' |sed 's/a/T/g'|sed 's/t/A/g'|sed 's/c/G/g'|sed 's/g/C/g' > a
paste <(cut -f 1 transcript.putativecds-.fa) <(less -S a|sed 's/ /\t/g'|cut -f 1)  > b
mv b transcript.putativecds-.fa
cat transcript.putativecds+.fa transcript.putativecds-.fa > transcript.putativecds.fa
perl /mnt/share/lisx/pacbio/bin/Base.faConvertProtein.fa.pl --Fasta <(less -S transcript.putativecds.fa |awk '{print ">"$0}'|sed 's/\t/\n/g' ) --output output.trans
### makedb for aa sequence of bgm transcriptome
mkdir blastdb_pro && cd blastdb_pro
cp ../output.trans transcript.fa
makeblastdb -in ../output.trans -dbtype prot -parse_seqids -out transcript.fa > blastdb.log 
### run blastp for refseq protein and genbank protein
blastp -query /mnt/share/lisx/pacbio/NCBI_inputdata/protein/refseq/human.refseq.aa.fa -out blastp.human.out -db blastdb_pro/transcript.fa -outfmt 7 -evalue 1e-5 -num_descriptions 10 -num_threads 10 > result/blastp.human.err 2>result/blastp.human.log &
blastp -query /mnt/share/lisx/pacbio/NCBI_inputdata/protein/refseq/rhesus.refseq.aa.fa -out blastp.rhesus.out -db blastdb_pro/transcript.fa -outfmt 7 -evalue 1e-5 -num_descriptions 10 -num_threads 10 > blastp.rhesus.err 2>blastp.rhesus.log &
blastp -query /mnt/share/lisx/pacbio/NCBI_inputdata/genbank/protein/bin-new/genbank.*.fa -out blastp.genebank.*.out -db ../../blastdb_pro/transcript.fa -outfmt 6 -evalue 1e-5 -num_threads 5 2>result/blastp.genebank.*.out.log &
