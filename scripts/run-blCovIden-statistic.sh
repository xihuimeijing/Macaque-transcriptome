#!/bin/sh
######   bash run-statistic.sh ${name} 
name=$1

###### gene && trans number:
less -S ${name} |cut -f 1|sort|uniq|awk '{a=a+1}END{print "GeneNumber:\t"a}' > statis.$name
less -S ${name} |cut -f 2|sort|uniq|awk '{a=a+1}END{print "TransNumber:\t"a}' >> statis.$name
###### pro-coverage of gene
echo -ne "\n\npro-coverage of gene:\n" >> statis.$name
paste <( perl /mnt/share/lisx/scripts/bin_value.pl -f <(perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f <(less -S ${name} |cut -f 1,3|sort|uniq)|sed 's/,/\t/g'|awk '{max=$2; if (NF>2){for (i=2;i<=NF;i++){if ($i>=max ){max=$i}} print $1"\t"max*10};if (NF==2){print $1"\t"$2*10} }' ) -value 2 -n 11 -max 11 -min 0 -equal min -outbin min~max -outnumber percent |awk '{a=a+$2;print $0"\t"a}'|sed 's/~/\t/g'|awk '{print $1/10"~"$2/10"\t"$3"\t"$4}'|sed 's/1~1.1/=1/g' ) <(perl /mnt/share/lisx/scripts/bin_value.pl -f <(perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f <(less -S ${name} |cut -f 1,3|sort|uniq)|sed 's/,/\t/g'|awk '{max=$2; if (NF>2){for (i=2;i<=NF;i++){if ($i>=max ){max=$i}} print $1"\t"max*10};if (NF==2){print $1"\t"$2*10} }' ) -value 2 -n 11 -max 11 -min 0 -equal min -outbin min~max -outnumber number|awk '{a=a+$2;print $0"\t"a}' |cut -f 2,3) >> statis.$name
####### pro-coverage of trans
echo -ne "\n\npro-coverage of trans:\n" >> statis.$name
paste <( perl /mnt/share/lisx/scripts/bin_value.pl -f <(perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f <(less -S ${name} |cut -f 2,3|sort|uniq)|sed 's/,/\t/g'|awk '{max=$2; if (NF>2){for (i=2;i<=NF;i++){if ($i>=max ){max=$i}} print $1"\t"max*10};if (NF==2){print $1"\t"$2*10} }' ) -value 2 -n 11 -max 11 -min 0 -equal min -outbin min~max -outnumber percent |awk '{a=a+$2;print $0"\t"a}'|sed 's/~/\t/g'|awk '{print $1/10"~"$2/10"\t"$3"\t"$4}'|sed 's/1~1.1/=1/g' ) <(perl /mnt/share/lisx/scripts/bin_value.pl -f <(perl /mnt/share/lisx/scripts/bin-mergeValue.pl -f <(less -S ${name} |cut -f 2,3|sort|uniq)|sed 's/,/\t/g'|awk '{max=$2; if (NF>2){for (i=2;i<=NF;i++){if ($i>=max ){max=$i}} print $1"\t"max*10};if (NF==2){print $1"\t"$2*10} }' ) -value 2 -n 11 -max 11 -min 0 -equal min -outbin min~max -outnumber number|awk '{a=a+$2;print $0"\t"a}' |cut -f 2,3) >> statis.$name

########## L1 L2 L3 Single
echo -ne "\n\nTerm\tall\tL1\tL2\tL3\tSingle\n" >> statis.$name
less -S ${name} |cut -f 1|sort|uniq|awk '{a=a+1}END{printf "Gene\t"a}' >> statis.$name
less -S ${name} |awk '{if ($2~/L1/ )print $0}' |cut -f 1|sort|uniq|awk '{a=a+1}END{printf "\t"a}' >> statis.$name
less -S ${name} |awk '{if ($2~/L2/ )print $0}'|cut -f 1|sort|uniq|awk '{a=a+1}END{printf "\t"a}' >> statis.$name
less -S ${name} |awk '{if ($2~/L3/ )print $0}'|cut -f 1|sort|uniq|awk '{a=a+1}END{printf "\t"a}' >> statis.$name
less -S ${name} |awk '{if ($2~/Single/ )print $0}'|cut -f 1|sort|uniq|awk '{a=a+1}END{printf "\t"a"\n"}' >> statis.$name

less -S ${name} |cut -f 2|sort|uniq|awk '{a=a+1}END{printf "Trans\t"a}' >> statis.$name
less -S ${name} |awk '{if ($2~/L1/ )print $0}'|cut -f 2|sort|uniq|awk '{a=a+1}END{printf "\t"a}' >> statis.$name
less -S ${name} |awk '{if ($2~/L2/ )print $0}'|cut -f 2|sort|uniq|awk '{a=a+1}END{printf "\t"a}' >> statis.$name
less -S ${name} |awk '{if ($2~/L3/ )print $0}'|cut -f 2|sort|uniq|awk '{a=a+1}END{printf "\t"a}' >> statis.$name
less -S ${name} |awk '{if ($2~/Single/ )print $0}'|cut -f 2|sort|uniq|awk '{a=a+1}END{printf "\t"a"\n"}' >> statis.$name

