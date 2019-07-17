#!/bin/sh

rm -r all-ctss-bed-peak && mkdir all-ctss-bed-peak 
echo -ne "SRR\tGEO\ttissue\n" > statistic.txt
sed -n '29p' SRR4281247/statisticReads.txt >> statistic.txt
echo -ne "ctss_number\n" > ctss.txt
for ((i=247;i<310;i++ ))
do
	grep SRR4281${i} name |sed 's/\n//g' >> statistic.txt
	sed -n '30p' SRR4281${i}/statisticReads.txt >> statistic.txt
	ctss=` gunzip -c SRR4281${i}/ctss.bed.gz |awk '{a+=1} END {print NR}' `
	echo -ne "$ctss\n" >> ctss.txt
done
awk '{if(0==NR%2)printf("%s\n",$0);else printf("%s\t",$0)}' statistic.txt > a
 paste a ctss.txt > statistic.txt && rm ctss.txt
cd all-ctss-bed-peak
for ((i=247;i<310;i++ ))
do	
	name=`grep SRR4281${i} ../name |sed 's/\t/-/g' `	
	cp ../SRR4281${i}/ctss.bed.gz ctss-${name}-bed.gz	
done
/mnt/share/lisx/tools/DIP/dpi1/identify_tss_peaks.sh -g /share/data/chr.size/rheMac8.size -i 'all-ctss.bed.gz/*.bed.gz' -o DIP-out 2> DIP.log
gunzip -c DIP-out/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed.gz |awk '{if ($6=="+"){print $1"\t"$2"\t"$2+1}else{ print $1"\t"$3-1"\t"$3 }}' > peak.position
