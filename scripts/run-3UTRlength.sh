#!/bin/sh
perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f <(less -S /rd1/user/liym/transcriptome/PA/brain/PA.onGene.tsv|cut -f 3,4,5|sed 's/,/\t/g'|awk '{for (i=2;i<=(NF+1)/2;i++){print $1"\t"$i"\t"$(i+NF/2) }  }'|awk '{if ($3>=2)print $1"\t"$2";"$3}' )|sed 's/,/\t/g'|sed 's/;/,/g' > gene.tmp
less -S gene.tmp |sed 's/,/\t/g'|awk '{sum=0;for (i=1;i<=(NF-1)/2;i++ ){ sum=sum+$(i*2+1);}; for (i=1;i<=(NF-1)/2;i++ ){ print $1"\t"$(i*2)"\t"$(i*2+1)"\t"$(i*2+1)/sum }  }' > gene.APA.weight
rm gene.tmp
less -S /rd1/user/liym/transcriptome/PA/brain/PA.onGene.tsv |cut -f 3,4 > gene.APA3utr
perl ../APA-newgpe.pl -f2 <(awk '{print 0"\t"$0}' ../rh8.longestUtr ) -f1 gene.APA3utr |cut -f 2-|awk '{print 0"\t"$0}' > newAPA.gpe
perl /mnt/share/lisx/scripts/gpe-utr-cds-length.pl -f1 newAPA.gpe |sed 's/:/\t/g'|cut -f 3,6-|awk '{if ($4=="+"){print $13"\t"$4"\t"$6"\t"$1}else{print $13"\t"$4"\t"$5"\t"$1 } }' > gene.utr.length
perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f <(perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f1 <(cut -f 1,3- gene.utr.length ) -f2 gene.APA.weight -n 2 -e1 1 -e2 1 |cut -f 1,2,3,7|awk '{print $1"\t"$2"\t"$3*$4 }' |cut -f 1,3 ) |sed 's/,/\t/g'|awk '{sum=0; for (i=2;i<=NF;i++){sum=sum+$i };print $1"\t"sum }' > gene.3utr
rm gene.APA3utr
rm gene.utr.length
rm newAPA.gpe


