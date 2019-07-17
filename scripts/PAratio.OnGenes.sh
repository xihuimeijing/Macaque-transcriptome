#!/bin/sh
usage(){
    cat << EOF
    Description: Find PA based on cluster and assign them to genes and calculate the PA usage/ratio (PA supported reads/all reads assigned to the gene) for each gene.
    Author: Yumei Li, 2018/12/06
    Revision: Yumei Li, 2018/12/29, Revise the processing of genes in different chrom/strand.
    Usage: perl $0 -g <inFile.gpe> -r <inFile.bed12+> -o <outDir> 
    Output file: Gene_id readCount RPM
    Options:
        -g  FILE        Gene structure file in GPE format without bin column.
        -r  FILE        Mapped iso-seq reads in Bed12/bed12+ format.
        -o  FILE        The directory name for all the output files.
        -h --help       Print this help information.
EOF
    exit 0
}

[ $1 ] || usage
outDir=$PWD
while getopts "hg:r:o:" OPTION
do
    case $OPTION in
        h) usage;;
        g) geneStructure=$OPTARG;;
        r) reads=$OPTARG;;
        o) outDir=$OPTARG;;
        ?) usage;;
    esac
done
perlScript='/rd1/user/liym/transcriptome/scripts/others'
perlMyself='/rd1/user/liym/transcriptome/scripts'
#1. Assign reads to genes
perl $perlScript/readsAssigner.pl -g $geneStructure -s $reads >$outDir/assign.Gene.bed12+
#2. Classify reads into unique, ambiguous or novel based on assigned gene number
perl $perlMyself/readGroup_YM.pl -g $geneStructure $outDir/assign.Gene.bed12+ >$outDir/unambi.bed12+ 2>$outDir/ambiguous.bed12+
#3. PA cluster for each gene
perl $perlMyself/paCluster_revision.pl -d 30 -m mode -w 3 $outDir/unambi.bed12+ >$outDir/paCluster.30size.bed8+
sort -k6,6 -k1,1 -k9,9 $outDir/paCluster.30size.bed8+ >$outDir/tmp;mv $outDir/tmp $outDir/paCluster.30size.bed8+
#4. Summary PA for each gene
firstGname=$(awk '$9!~","' $outDir/paCluster.30size.bed8+ | head -n1 |cut -f9)
awk '$9!~","' $outDir/paCluster.30size.bed8+ | awk -v OFS="\t" -v gene=$firstGname '
        BEGIN{geneName=gene;PA="";readNum="";sum=0;chr="";strand=""}
        {
            if($9 == geneName){
                if(sum==0){
                    geneName=$9;chr=$1;strand=$6;sum+=$5;PA=$8;readNum=$5;
                }else{
                    sum+=$5;PA=PA","$8;readNum=readNum","$5;
                }
            }else{
                if(readNum~","){
                    split(readNum,a,",");ratio=sprintf("%.4f",a[1]/sum);for(i=2;i<=length(a);i++){ratio=ratio","sprintf("%.4f",a[i]/sum)}
                }else{
                    ratio=1.00 
                }
                print chr,strand,geneName,PA,readNum,ratio;
                geneName=$9;chr=$1;strand=$6;sum=$5;PA=$8;readNum=$5;
            }
        }
        END{
            if(readNum~","){
                split(readNum,a,",");ratio=sprintf("%.4f",a[1]/sum);for(i=2;i<=length(a);i++){ratio=ratio","sprintf("%.4f",a[i]/sum)}
            }else{
                ratio=1.00 
            }
            print chr,strand,geneName,PA,readNum,ratio;
        }' >$outDir/PA.onGene.tsv
