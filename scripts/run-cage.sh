#!/bin/sh
SQgpe=/mnt/share/shenq/forSX/20180517/brain.gpe
tissue=brain
paste <(cut -f 1 $SQgpe |sed 's/-/_/g' ) <(cut -f 2- $SQgpe) |sed 's/brain_//g'|sed 's/brain-//g' > ${tissue}.cds.gpe
bash runSR-.sh ${tissue}
bash runSR+.sh ${tissue}
cat *level1* *level2* > trans.gene.SRcage
rm *level*
bash runGR-.sh ${tissue}
bash runGR+.sh ${tissue}
cat *level3* *level4* > trans.gene.GRcage
rm *level*
perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f1 trans.gene.SRcage -f2 trans.gene.GRcage -n 1 -e1 3 -e2 3 -only only f2 > trans.gene.GRcage.tmp 
cat trans.gene.GRcage.tmp trans.gene.SRcage  > trans.gene.levelCage
mkdir test
mv trans.gene.*cage* test/

############# gpe file of different level
perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f1 trans.gene.levelCage -f2 brain.cds.gpe -n 1 -e1 3 -e2 17 |grep -v level2|grep -v level4 |cut -f 1,3,6-|awk '{print $1"."$2"\t"$0}'|cut -f 1,4- > brain-level.tmp1,3.gpe

perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f1 trans.gene.levelCage -f2 brain.cds.gpe -n 1 -e1 3 -e2 17 |grep -v level1|grep -v level3 |cut -f 1,3,4,6-|awk '{print $1"."$2"\t"$0}'|cut -f 1,4-|grep + > brain-level.tmp24.gpe+
less -S brain-level.tmp24.gpe+ |cut -f 2 |sed 's/-/\t/g'|sed 's/:/\t/g'|cut -f 3 > brain-level.tmp24.gpe+.start
paste brain-level.tmp24.gpe+.start <(less -S brain-level.tmp24.gpe+|cut -f 10|sed 's/,/\t/g'|cut -f 2- ) |sed 's/\t/,/g' > brain-level.tmp24.gpe+.startexon
paste brain-level.tmp24.gpe+.start brain-level.tmp24.gpe+.startexon brain-level.tmp24.gpe+|awk '{$7=$1;$12=$2;print $0}'|sed 's/ /\t/g'|cut -f 3,5- > brain-level.tmp2,4.gpe+
rm brain-level.tmp24.gpe+*

perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f1 trans.gene.levelCage -f2 brain.cds.gpe -n 1 -e1 3 -e2 17 |grep -v level1|grep -v level3 |cut -f 1,3,4,6-|awk '{print $1"."$2"\t"$0}'|cut -f 1,4-|grep -v + > brain-level.tmp24.gpe-
less -S brain-level.tmp24.gpe- |cut -f 2 |sed 's/:/\t/g'|cut -f 3|sed 's/-/\t/g'|sed 's/,/\t/g'|cut -f 2 > brain-level.tmp24.gpe-.start
paste <(less -S brain-level.tmp24.gpe-|cut -f 11|sed 's/,/\t/g'|awk '{for (i=1;i<NF;i++){printf $i"\t"}; printf "\n"}') brain-level.tmp24.gpe-.start |sed 's/\t\t/\t/g'|sed 's/^\t//g'|sed 's/\t/,/g'|awk '{print $0","}' > brain-level.tmp24.gpe-.startexon
paste brain-level.tmp24.gpe-.start brain-level.tmp24.gpe-.startexon brain-level.tmp24.gpe-|awk '{$8=$1;$9=$1;$10=$1;$13=$2;print $0}'|sed 's/ /\t/g'|cut -f 3,5- > brain-level.tmp2,4.gpe-
rm brain-level.tmp24.gpe-*
cat brain-level.tmp1,3.gpe brain-level.tmp2,4.gpe*|sed 's/,,/,/g' > brain-level.gpe
rm *tmp*
gpe2bed.pl brain-level.gpe > brain-level.bed

######### brain cds+cage level gpe file
cat brain-level.gpe <(perl /mnt/share/lisx/scripts/comm_file1_file2.pl -f1 brain.cds.gpe -f2 <(less -S brain-level.gpe |sed 's/\./\t/g' ) -n 1 -e1 16 -e2 16 -only onlyf1 ) |awk '{print "brain_"$0}' > brain-cds-level.gpe
paste <(cut -f 1-11 brain-cds-level.gpe ) <( awk '{print "brain-"$12}' brain-cds-level.gpe ) <(cut -f 13- brain-cds-level.gpe ) > a
mv a brain-cds-level.gpe
awk '{print "heart_"$1"\theart-"$2"\t"$3"\t"$4}'  trans.gene.levelCage > a
mv a trans.gene.levelCage
