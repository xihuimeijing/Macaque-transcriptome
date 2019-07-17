#!/bin/sh
#for (( i=259;i<311;i++ ))
for ((i=247;i<310;i++ ))
do
#i=249
	wget ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR428/SRR4281${i}/SRR4281${i}.sra
	mkdir SRR4281${i}  
	cd SRR4281${i}
	fastq-dump ../SRR4281${i}.sra --split-files -O .
	mkdir fastqc_out && fastqc -o fastqc_out -f fastq --extract SRR4281${i}_1.fastq SRR4281${i}_2.fastq 2> fastqc_out/fastqc.log 
	rm ../SRR4281${i}.sra
	less -S fastqc_out/SRR4281${i}_1_fastqc/fastqc_data.txt |grep 'Index'|cut -f 1 > index5
	less -S fastqc_out/SRR4281${i}_2_fastqc/fastqc_data.txt |grep 'Primer'|cut -f 1 > primer3
	cp SRR4281${i}_1.fastq precutadapte-SRR4281${i}_1.fastq 
	cp SRR4281${i}_2.fastq precutadapte-SRR4281${i}_2.fastq
	for ((ii=1;ii<=2;ii++ ))
	do
		Index5=`sed -n "$ii"'p' index5 `	
		Primer3=` sed -n "$ii"'p' primer3 `
		cutadapt -a $Index5 -A $Primer3 -o SRR4281${i}_1_cutadapter.fastq -p SRR4281${i}_2_cutadapter.fastq -m 20 precutadapte-SRR4281${i}_1.fastq precutadapte-SRR4281${i}_2.fastq > cutAdapter${ii}.log 
		mv SRR4281${i}_1_cutadapter.fastq precutadapte-SRR4281${i}_1.fastq && mv SRR4281${i}_2_cutadapter.fastq precutadapte-SRR4281${i}_2.fastq
	done	
	mv precutadapte-SRR4281${i}_1.fastq SRR4281${i}_1_cutadapter.fastq && mv precutadapte-SRR4281${i}_2.fastq SRR4281${i}_2_cutadapter.fastq
	fqTrimer.pl -l [ATCGN]{6} SRR4281${i}_1_cutadapter.fastq >SRR4281${i}_1_trim6base.fastq 
	fqTrimer.pl -r [ATCGN]{9} SRR4281${i}_2_cutadapter.fastq >SRR4281${i}_2_trim6base.fastq
	fqPeFilter.pl -f 0.02 -b N -q 0.2 -l 25 -a 25 -t 5,20 -o1 SRR4281${i}_fin_1.fq -o2 SRR4281${i}_fin_2.fq  SRR4281${i}_1_trim6base.fastq SRR4281${i}_2_trim6base.fastq > qulity_filter.log
	perl /mnt/share/shenq/tool.pl/05_Capseq/CapSeqReadFilter.pl -r1 SRR4281${i}_fin_1.fq -r2 SRR4281${i}_fin_2.fq -f1 3 -f2 0 -o1 SRR4281${i}_final_1.fq -o2 SRR4281${i}_final_2.fq -minlength 20 -removelist no3Greads 
	mkdir fastqc_out_final && fastqc -o fastqc_out_final -f fastq --extract SRR4281${i}_final_1.fq SRR4281${i}_final_2.fq 2> fastqc_out_final/fastqc.log 
	rm *.fastq
	rm SRR4281${i}_fin_*.fq
	tophat2 --read-mismatches 8 --read-gap-length 3 --read-edit-dist 10 --min-anchor 8 --num-threads 35 --mate-inner-dist 0 --no-coverage-search --segment-mismatches 2 --library-type fr-unstranded /share/data/bowtie2/rheMac8/rheMac8 SRR4281${i}_final_1.fq SRR4281${i}_final_2.fq > tophat2.log 2> tophat2.err 
	/home/share/local/bin/samtools stats tophat_out/accepted_hits.bam > bamstats_accepted.txt 
################### statistic 
	
	read1N=` less -S fastqc_out/SRR4281${i}_1_fastqc/fastqc_data.txt |grep 'Total Sequences' |sed 's/Total Sequences//g'|sed 's/^[ \t]*//g' `
	read2N=` less -S fastqc_out/SRR4281${i}_2_fastqc/fastqc_data.txt |grep 'Total Sequences' |sed 's/Total Sequences//g'|sed 's/^[ \t]*//g' `	
	let "Allread=read1N+read2N"

	read1n=` less -S fastqc_out_final/SRR4281${i}_final_1_fastqc/fastqc_data.txt |grep 'Total Sequences' |sed 's/Total Sequences//g'|sed 's/^[ \t]*//g' `
	read2n=` less -S fastqc_out_final/SRR4281${i}_final_2_fastqc/fastqc_data.txt |grep 'Total Sequences' |sed 's/Total Sequences//g'|sed 's/^[ \t]*//g' `
	let "allread=read1n+read2n"	
	let "read1f=read1N-read1n"
	let "read2f=read2N-read2n"
	echo -ne "raw reads counts\nsequence1\t${read1N}\tfiltered\t${read1f}\nsequence2\t$read2N\tfiltered\t${read2f}\nallsequence\t$Allread\n" > statisticReads.txt
	echo -ne "\nprepare for mapping reads counts\nsequence1\t$read1n\nsequence2\t$read2n\nallsequence\t$allread\n" >> statisticReads.txt	
####################	
	samtools view -bu -q 50 -@ 5 tophat_out/accepted_hits.bam | samtools sort -@ 5 -m 10G - uniq.sorted 2> samtoolsSort.log  
	mv uniq.sorted.bam tophat_out/  
	/home/share/local/bin/samtools stats tophat_out/uniq.sorted.bam > bamstats_uniq.txt 
	samtools view -bu -f 65 -F 8 -@ 5 tophat_out/accepted_hits.bam | samtools sort -@ 5 -m 10G - tophat_out/read1.sorted 
	/home/share/local/bin/samtools stats tophat_out/read1.sorted.bam > bamstats_read1.txt

	echo -ne "\n\n######### mapping result ###############\n" >> statisticReads.txt
	cat tophat_out/align_summary.txt >> statisticReads.txt 	
	readmapping=` grep "raw total sequences" bamstats_accepted.txt |sed 's/raw total sequences://g'|sed 's/SN//g'|sed 's/^[ \t]*//g' `
	readuniq=` grep "raw total sequences" bamstats_uniq.txt |sed 's/raw total sequences://g'|sed 's/SN//g'|sed 's/^[ \t]*//g' `
	readsequence1=` grep "raw total sequences" bamstats_read1.txt |sed 's/raw total sequences://g'|sed 's/SN//g'|sed 's/^[ \t]*//g' `
	filter_percent=`awk 'BEGIN{printf '$read1f'/'$read1N'}'`
	let filter=read1f*2
	mapping_percent=`awk 'BEGIN{printf '$readmapping'/'$allread'}'`
	uniq_percent=` awk 'BEGIN{printf '$readuniq'/'$readmapping'}' `
	echo -ne "\n############ summary ###########\nall_reads\tfitter_reads\tfilter_percent\treads_PreMap\taccepted_reads\tmapping_percent\tuniq_reads\tuniq_percent\tread1\n${Allread}\t${filter}\t${filter_percent}\t${allread}\t${readmapping}\t${mapping_percent}\t${readuniq}\t${uniq_percent}\t${readsequence1}\n" >> statisticReads.txt
	bedtools bamtobed -bed12 -i tophat_out/read1.sorted.bam > read1.bed12 
	less -S read1.bed12 |awk '{if ($6=="+"){print $1"\t"$2"\t"$6}else{print $1"\t"$3-1"\t"$6} }'|awk '{print $1"#"$2"#"$2+1"#"$1":"$2".."$2+1","$3"#"$3}'|sort|uniq -c |sed 's/^[ \t]*//g'|awk '{print $2"\t"$1}' |sed 's/#/\t/g'|awk '{print $0"\t"$5}'|cut -f 1-4,6- |sort -k 1,1 -k 2,2n | gzip - > ctss.bed.gz
	cd ../
done

