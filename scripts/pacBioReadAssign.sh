#!/usr/bin/env bash

usage(){
  cat << EOF
Description:
    Output the pacbio reads number for each input reference transcript.
Usage:
    sh $0 -g gene.bed12 -r reads.bed12 -p tmp_prefix -o rst.bed12+
Output: bed12_fields_reads transcritpName
Author: Yumei Li, 2018/07/27
Options:
    -g  FILE          The gene structure file in BED12 format.(If not provided, can be read from STDIN)
    -r 	FILE          Mapped pacBio reads in BED12 format.
    -p  STRING        The prefix for the intermediate file.
    -o  FILE          The output file name.
    -h --help         Print this help information
EOF
    exit 0
}

[ $1 ] || usage

while getopts "hg:r:p:o:" OPTION
do
    case $OPTION in
        h) usage;;
        g) geneStructure=$OPTARG;;
        r) readsFile=$OPTARG;;
        p) prefix=$OPTARG;;
        o) outFile=$OPTARG;;
        ?) usage;;
    esac
done
shift $((OPTIND - 1))

if [ -z $geneStructure ];then
  bedtools intersect -loj -s -a stdin -b $readsFile >${prefix}.intersect.bed
else
  bedtools intersect -loj -s -a $geneStructure -b $readsFile >${prefix}.intersect.bed
fi
awk '$14>=0' ${prefix}.intersect.bed|cut -f1-12|bedToGenePred stdin ${prefix}.trans.gpe
awk '$14>=0' ${prefix}.intersect.bed|cut -f13-24|bedToGenePred stdin ${prefix}.read.gpe
overlapLine=$(wc -l ${prefix}.trans.gpe|awk '{print $1}')
if [ $overlapLine -gt 0 ];then
  paste ${prefix}.trans.gpe ${prefix}.read.gpe|awk -v OFS="" -v ORS="" '
    {
      if($18>1){
        split($19,s,",");
        startR=s[2];
        split($20,e,",");
        endR=e[1];
        for(j=3;j<=$18;j++){startR=startR","s[j]}
        for(k=2;k<=$18-1;k++){endR=endR","e[k]}
        if($9~startR && $10~endR){
          for(i=11;i<=20;i++){
            print $i"\t";
          }
          print $1"\n";
        }
      }
    }' >$outFile
  rm ${prefix}.intersect.bed ${prefix}.trans.gpe ${prefix}.read.gpe
else
  touch $outFile
fi