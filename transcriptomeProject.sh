#!/bin/sh
#All commands for transcriptome project

#1 Iso-seq data processing (@venus)
	##PATH set: you should have python2.7 biopython2.7 phmmer in your corresponding PATH  
	set -x 
	##parameters
	list=$1
	thread=$2
	minCCSlength=50
	species=rheMac8
	tissue=brain
	npass=5
	##scripts
	smrt_dir=/rd1/user/shenq/smrtlink/smrtlink/smrtcmds/bin
	cDNA_dir=/rd1/user/shenq/Annotation/bin/cDNA_primer-master/scripts
	script_dir=/rd1/user/shenq/Annotation/bin
	cdsLength=300 #100 aa
	##files
	primer=/rd1/user/shenq/smrtlink/smrtlink/mydata/primers.fa #
	ref=/rd1/user/shenq/reference/fna/rheMac8.fa
	gmap=/rd1/user/shenq/reference/gmap/rheMac8/rheMac8
	gap=/rd1/user/shenq/reference/gap/rheMac8.bed4+
	twoBit=/rd1/user/shenq/reference/2Bit/rheMac8_softMasked.2bit
	rnaSeqJunc=/rd1/user/shenq/Annotation/RNAseq/juncFile/RBrhesus/All.RNAseq.MT3.junc.uniq
	rnaSeqIntron=/rd1/user/shenq/Annotation/RNAseq/juncFile/RBrhesus/RNAseq.MT3.intron
	rnaSeqIntronBed12=/rd1/user/shenq/Annotation/RNAseq/juncFile/RBrhesus/RNAseq.sam2junc.bed12.MT3
	
	##1.1 Extract CCS(ROI)
		###dataset merge
		file=`perl -ne 's/\n/ /g;print;' $list`
		$smrt_dir/dataset merge file.subreadset.xml $file
		###ccs extraction
		$smrt_dir/ccs --polish --maxLength 30000 --minLength 50 --minPasses 0 --minPredictedAccuracy 0.75 --minZScore -9999 --maxDropFraction 0.8 --minSnr 3.75 --minReadScore 0.75 --numThreads $thread file.subreadset.xml  ccs.bam 2>>ccs.log
		$smrt_dir/dataset create --type ConsensusReadSet ccs.consensusreadset.xml ccs.bam
		$smrt_dir/bam2fastq ccs.bam -u -o ccs
	
	##1.2 Summary quality
		samtools view ccs.bam |perl -ne '/\tnp:i:(.*)\trq:f:(.*)\trs/;@a=split/\t/;$b=length($a[9]);print "$a[0]\t$1\t$b\t$2\n";' >npass.ccs.txt
		cut -f 2 npass.ccs.txt|$script_dir/hist.R -x='Pass Number' -y=Density -b=1 -d -x1=0 -x2=50 -p=PassNumber.pdf 2>/dev/null
		cut -f 4 npass.ccs.txt|$script_dir/hist.R -x='CCS Quality' -y=Density -b=100 -d -x1=0 -x2=1 -p=ccsQuality.pdf 2>/dev/null
		awk '$2>="$npass"' npass.ccs.txt|cut -f 3 |$script_dir/hist.R -x='Reads Length (greater than npass)' -y=Density -b=100 -d -x1=0 -x2=10000 -p=PassNumberGreaterLength.pdf 2>/dev/null
	
	##1.3 Get non chimeric FL reads
		$cDNA_dir/hmmer_wrapper.py --primer_search_window 150 --min-score 8 --min-seqlen $minCCSlength --cpus $thread --change-seqid --primer_filename $primer --input_filename ccs.fastq --output-fq --output_filename FL.fq --directory hmmer/ --must-see-polyA --left-nosee-ok 2> hmmer_wrapper.log
		$script_dir/fqFormater.pl -t ccs.fastq |perl -ne 's/\/ccs//;print;'|$script_dir/filter.pl -o <($script_dir/fqFormater.pl -t FL.fq | cut -f1 | sed 's/\/[0-9_]*CCS$//') |$script_dir/fqFormater.pl -f >noFL.fq 
		
		$cDNA_dir/chimera_finder.py --min_dist-from_end 50 --cpus $thread --primer_filename $primer --input_filename FL.fq --directory chimera_FL/ --output-fq 2>chimera_finder.FL.log \
		&& mv FL.fq.non_chimera.fq FLnoC.fq && mv FL.fq.is_chimera.fq FLC.fq
		$cDNA_dir/chimera_finder.py --min_dist-from_end 50 --cpus $thread --primer_filename $primer --input_filename noFL.fq --directory chimera_noFL/ --output-fq 2>chimera_finder.noFL.log \
		&& mv noFL.fq.non_chimera.fq noFLnoC.fq && mv noFL.fq.is_chimera.fq noFLC.fq
	
	##1.4 Trim PA tail
		mkdir -p recovery
		for f in $(seq 0.4 0.05 0.7);do for w in $(seq 30 1 60);do recoveryLine=$($script_dir/separatePaReadsFromNotFL.pl noFLnoC.fq -f $f -w $w 2>/dev/null | wc -l); recoveryReads=$(($recoveryLine/4)); echo -e "$w\t$recoveryReads" >>recovery/f=$f.lst;done ;done && $script_dir/lines.R -x='Window Size' -y=Recovery -legendT=Fraction -w=15 recovery/*.lst -p=recovery.pdf 2>/dev/null
		w=40
		f=0.65
		$script_dir/separatePaReadsFromNotFL.pl noFLnoC.fq -w $w -f $f >recovery.fq 2>separatePaReadsFromNotFL.log
		($script_dir/fqTrimPA.pl -w $w -f $f FLnoC.fq 2>fqTrimPA.log; cat recovery.fq) | $script_dir/fqSeFilter.pl --minLength $minCCSlength >paTrimmed.fq
		perl $script_dir/fqLength.pl paTrimmed.fq >paTrimmed.Length
		cut -f 2 paTrimmed.Length|$script_dir/hist.R -x='Reads Length (paTrimmed)' -y=Density -b=100 -d -x1=0 -x2=10000 -p=PaTrimmed.ReadLength.pdf 2>/dev/null
		ln -s paTrimmed.fq final.fq
		###PA length summary
		(awk '$3<$2' fqTrimPA.log | sed 's/^@//' | $script_dir/skyjoin - FL.fq.primer_info.txt | awk '$10>$9{print $1"\t"$10-$9+$2-$3}'
		awk '$3==$2' fqTrimPA.log | sed 's/^@//' | $script_dir/filter.pl -o - -m i FL.fq.primer_info.txt | awk '$7!="NA" && $8!="NA" && $8-$7>0{print $1"\t"$8-$7}'
		awk '$3<$2{print $1"\t"$2-$3}' separatePaReadsFromNotFL.log| sed 's/^@//' ) | sort -n >paTailLength.tsv
		cut -f2 paTailLength.tsv | $script_dir/hist.R -x='Tail length' -y=Density -b=1 -d -x1=0 -x2=100 -p=paTailLength.pdf 2>/dev/null
	
	##1.5 Mapping
		mkdir -p gmap && cd gmap
		$smrt_dir/gmap -D $gmap -d gmap_db --input-buffer-size=10000 --batch=5 --min-intronlength=70 -K 1100000 --totallength=2500000 --nthreads=$thread --format=samse --ordered --split-output=sam --output-buffer-size=10000 --read-group-name=$species --read-group-library=$tissue --read-group-platform=PacBio ../final.fq 2> gmap.log
		(grep "^@" sam.uniq; grep -vh "^@" sam.*) |$script_dir/samAddTag.pl --coverage --identity --unmapped unmapped.sam 2>lengthInconsistent.sam | samtools view -buS - |samtools sort -m 15G -@ $thread -o mapped.sorted.bam - && $script_dir/sam2bed.pl -s -t CV,ID mapped.sorted.bam >mapped.bed12+
		$script_dir/readsFilter.pl -r 0.8 mapped.bed12+ 2>non-uniq.bed12+ >uniq.bed12+
		(samtools view -H mapped.sorted.bam; samtools view mapped.sorted.bam |$script_dir/filter.pl -o <(awk 'BEGIN{FS=OFS="\t"}{print $1,$2+1,$4,$5}' uniq.bed12+) -1 1,2,3,4 -2 3,4,1,5 -m i) | samtools view -buS - | samtools sort -m 8G -@ $thread -o uniq.sorted.bam -
		$script_dir/bed2gpe.pl uniq.bed12+ | $script_dir/gpeFeature.pl -e | bedtools intersect -a - -b $gap -wa -u | $script_dir/filter.pl -o /dev/stdin -1 7 -2 4 uniq.bed12+ >noGap.bed12+
		$script_dir/dnaOrInternalPrimingContFilter.pl -b $twoBit -t <(sed 's/^@//' ../paTailLength.tsv) -r dnaOrInternalPrimingContFilter.log noGap.bed12+ >deCont.bed12+ 2>dnaCont.bed12+
		(samtools view -F 0x10 uniq.sorted.bam | perl -ne 'print if (split "\t")[5] =~ /(\d+)S$/ && $1 >30'; samtools view -f 0x10 uniq.sorted.bam | perl -ne 'print if (split "\t")[5] =~ /^(\d+)S/ && $1 >30' ) | $script_dir/filter.pl -o /dev/stdin -2 4 deCont.bed12+ >processed.bed12+
		###Evaluation of error rate
		awk 'BEGIN{OFS="\t"}{print $4,100-$13*$14/100}' processed.bed12+ | tee errorRate.log | cut -f2 |$script_dir/distrCurve.R -d -m='Rate of Mismatch & Indel for Mapping' -x=Rate -y='Fraction of Reads' -p=errorRate.pdf 2>/dev/null
		###Format the results
		cat processed.bed12+ |perl -ne  's/\/[0-9_]*CCS/\/ccs/;print;' >processed.bed12+.format
		perl $script_dir/select.pl -i <(awk -v npass="$npass" '$2>=npass' ../npass.ccs.txt) -r processed.bed12+.format -c1 1 -c2 4 >processed.npass.bed12
		perl $script_dir/select.pl -i <(awk '$2>=1' ../npass.ccs.txt) -r processed.bed12+.format -c1 1 -c2 4 >processed.pass1.bed12
		cd ..
		cat gmap/uniq.bed12+ |bedToGenePred stdin stdout|perl $script_dir/length.pl -i - |cut -f 2|$script_dir/hist.R -x='Reads Length (Uniq)' -y=Density -b=100 -d -x1=0 -x2=10000 -p=uniqLength.pdf 2>/dev/null &
#2 The definition of gene models
	##2.1 Get multi- and single-exon isoforms
		mkdir -p isoform && cd isoform
		ln -s ../gmap/processed.npass.bed12 && ln -s ../gmap/processed.pass1.bed12 && ln -s ../gmap/uniq.bed12+
		awk '$10>1' processed.npass.bed12|bedToGenePred stdin processed.npass.multiExon.gp
		bedToGenePred  processed.pass1.bed12 processed.pass1.gp
		awk '$8==1' processed.pass1.gp > processed.pass1.singleExon.gp
		awk '$8>1' processed.pass1.gp|$script_dir/filter.pl -o processed.npass.multiExon.gp -1 1 -2 1 - > processed.pass1npass.multiExon.gp
	##2.2 Define isoforms in different levels
		###Multi-exon isoforms filtered by RNA-seq junctions (at least 2 junc-reads supported junction)
		perl -alne 'if ($F[7]>1) {@a=split/,/,$F[8];@b=split/,/,$F[9];shift @a;pop @b;$F[8]=join ",",@a;$F[9]=join ",",@b;print join "\t",@F;}' processed.npass.multiExon.gp > processed.npass.MultiexonMiddlejunc.tmp
		cut -f 2,3,8,9,10 processed.npass.MultiexonMiddlejunc.tmp|sort |uniq -c|perl -ne 's/(  )( )*//g;s/ /\t/g;print ' >processed.npass.MultiexonMiddlejunc.count
		cat processed.npass.MultiexonMiddlejunc.count|$script_dir/hist.R -x='Identical Junction Count' -y=Density -b=1 -d -x1=0 -x2=25 -m=IdenticalJuncCounts -p=IdenticalJuncCounts.pdf 2>/dev/null
		$script_dir/filter.pl -o <(cat processed.npass.MultiexonMiddlejunc.tmp |cut -f 2,3,8,9,10|sort|uniq -d) -1 1,2,3,4,5 -2 2,3,8,9,10 -m i processed.npass.MultiexonMiddlejunc.tmp|cut -f 1 |perl $script_dir/select.pl -i - -r processed.npass.multiExon.gp -c1 1 -c2 1 >confirm.r1.gp.tmp1
		perl $script_dir/select.pl -i confirm.r1.gp.tmp1 -r processed.npass.multiExon.gp -c1 1 -c2 1 -o|awk '$8>1' >non_confirm.r1.gp.tmp
		perl $script_dir/bin/gpe2Junction.pl -nb -i confirm.r1.gp.tmp1 |cut -f 1-3,5 |sort -u |perl -alne '$n++;$F[5]=$F[3];$F[4]=0;$F[3]="ConfirmR1Junc_"."$n";print join "\t",@F;' >confirm.r1.junc
		perl $script_dir/juncRNAseqConfirm.pl -j confirm.r1.junc -i non_confirm.r1.gp.tmp >confirm.r1.gp.tmp2
		perl $script_dir/select.pl -i confirm.r1.gp.tmp2 -r non_confirm.r1.gp.tmp -c1 1 -c2 1 -o >non_confirm.r1.gp
		cat confirm.r1.gp.tmp1 confirm.r1.gp.tmp2 >confirm.r1.gp
		##Rank
		perl $script_dir/juncRNAseqConfirm.pl -j <(grep newMonkey_R-"$tissue" $rnaSeqJunc) -i <(awk '$8>1' processed.pass1.gp) >confirm.r2.L.gp
		perl $script_dir/select.pl -i confirm.r2.L.gp -r confirm.r1.gp -c1 1 -c2 1 -o >confirm.r3.L.gp
		perl $script_dir/select.pl -i <(cat confirm.r2.L.gp confirm.r3.L.gp) -r <(awk '$8>1' processed.pass1.gp) -c1 1 -c2 1 -o >non_confirm.r2r3.L.gp
		genePredToBed non_confirm.r2r3.L.gp non_confirm.r2r3.L.bed
		##Use junc supported by R2 and RNA-seq junc to amend other isoform. In fact, just those RNA-seq junc make sense.
		perl $script_dir/removeMisAlignExon.pl -g confirm.r2.L.gp -j <(grep newMonkey_R-"$tissue" $rnaSeqIntronBed12) non_confirm.r2r3.L.bed | perl $script_dir/fillMissExon.pl -g confirm.r2.L.gp -j <(grep newMonkey_R-"$tissue" $rnaSeqIntronBed12) >realign.bed
		perl $script_dir/juncConsensus.pl -s <($script_dir/juncScoring_sq.pl -j <(grep newMonkey_R-"$tissue" $rnaSeqIntron) processed.pass1.bed12) realign.bed |bedToGenePred stdin amend.r3.gp.tmp #time consuming
		cat amend.r3.gp.tmp |$script_dir/gpeFeature.pl -e |awk '$8<0' |cut -f 7|$script_dir/filter.pl -o - -1 1 -2 1 amend.r3.gp.tmp > amend.r3.gp #remove exon-length <0 ;bug of above scripts
		rm *.tmp*
	##2.3 Remove redundant isoforms
		mkdir -p cuffcompare.all && cd cuffcompare.all
		awk '$10>1' ../uniq.bed12+ |cut -f 1-12 >ccs_uniqMapping.bed
		ln -s ../confirm.r2.L.gp && ln -s ../confirm.r3.L.gp && ln -s ../amend.r3.gp && ln -s ../processed.pass1.singleExon.gp
		cat confirm.r2.L.gp confirm.r3.L.gp amend.r3.gp processed.pass1.singleExon.gp >final.gp
		cat final.gp |genePredToGtf file stdin final.gtf
		cuffcompare -F final.gtf
		gtfToGenePred -genePredExt cuffcmp.combined.gtf cuffcmp.combined.gpe
		###level the transcript
		perl -alne 'if ($F[7]>1) {@a=split/,/,$F[8];@b=split/,/,$F[9];shift @a;pop @b;$F[8]=join ",",@a;$F[9]=join ",",@b;print join "\t",@F;}' cuffcmp.combined.gpe > all.tmp
		perl -alne 'if ($F[7]>1) {@a=split/,/,$F[8];@b=split/,/,$F[9];shift @a;pop @b;$F[8]=join ",",@a;$F[9]=join ",",@b;print join "\t",@F;}' confirm.r2.L.gp |cut -f 2,3,8,9,10|sort -u > L1.tmp
		perl -alne 'if ($F[7]>1) {@a=split/,/,$F[8];@b=split/,/,$F[9];shift @a;pop @b;$F[8]=join ",",@a;$F[9]=join ",",@b;print join "\t",@F;}' confirm.r3.L.gp |cut -f 2,3,8,9,10|sort -u> L2.tmp
		perl $script_dir/filter.pl -o L1.tmp -1 1-5 -2 2,3,8,9,10 -m i all.tmp |cut -f 1|perl $script_dir/select.pl -i - -r cuffcmp.combined.gpe -c1 1 -c2 1 |perl -lane '$F[0]=$F[0]."-L1";print join "\t",@F;' > level.gpe #iso-seq support
		perl $script_dir/filter.pl -o L1.tmp -1 1-5 -2 2,3,8,9,10 all.tmp >level.1.tmp
		perl $script_dir/filter.pl -o L2.tmp -1 1-5 -2 2,3,8,9,10 -m i level.1.tmp |cut -f 1|perl $script_dir/select.pl -i - -r cuffcmp.combined.gpe -c1 1 -c2 1 |perl -lane '$F[0]=$F[0]."-L2";print join "\t",@F;'  >> level.gpe #RNA-seq support
		perl $script_dir/filter.pl -o L2.tmp -1 1-5 -2 2,3,8,9,10 level.1.tmp |cut -f 1|perl $script_dir/select.pl -i - -r cuffcmp.combined.gpe -c1 1 -c2 1 |perl -lane '$F[0]=$F[0]."-L3";print join "\t",@F;' >>level.gpe #RNA-seq and Iso-seq R1 amend
		awk '$8==1' cuffcmp.combined.gpe|perl -lane '$F[0]=$F[0]."-Single";print join "\t",@F;' >>level.gpe #single exon
		###TooLong (not in L1 now)
		awk '$8>1' level.gpe>level.m.gpe
		perl $script_dir/length.pl -i level.m.gpe|cut -f 1,3 |join -1 1 -2 1 level.m.gpe - |perl -ne 's/ /\t/g;print;' |perl $script_dir/skycut.pl -f 2,4,5,1,8,3,12,16|perl -lane '$F[8]=$F[2]-$F[1];print join "\t",@F;' >level.m.bed.tmp
		bedtools intersect -s -wo -a level.m.bed.tmp -b level.m.bed.tmp| awk '$7!=$16'|awk '$18*0.9<=$19' |cut -f 1-9|sort -u |grep -v L1 >level.toolong.bed #Not in L1 now
		cat level.toolong.bed|cut -f 4 |cut -d '-' -f 1 |perl $script_dir/select.pl -i - -r cuffcmp.tracking -c1 1 -c2 1 | 
		cut -f 5 |cut -d '|' -f 1|cut -d '/' -f 1-2|cut -d ':' -f 2 |while read f;do grep "$f/" ../../gmap/processed.bed12+ ;done >level.toolong.reads.bed12+
		cut -f 14 level.toolong.reads.bed12+|$script_dir/hist.R -x='Identity' -y=Density -b=1 -d -x1=70 -x2=100 -m=toolong -p=TooLong_identity.pdf
		perl $script_dir/select.pl -i level.toolong.reads.bed12+ -r ../../gmap/processed.bed12+ -c1 4 -c2 4 -o |cut -f 14 |$script_dir/hist.R -x='Identity' -y=Density -b=1 -d -x1=70 -x2=100 -m=Nottoolong -p=NotTooLong_identity.pdf
		perl $script_dir/select.pl -i level.toolong.bed -r level.m.gpe -c1 4 -c2 1 |perl -lane '$F[0]=$F[0]."TooLong";print join "\t",@F' >level.gpe.tmp
		perl $script_dir/select.pl -i level.toolong.bed -r level.m.gpe -c1 4 -c2 1 -o >>level.gpe.tmp
		###Single
		awk '$8==1' level.gpe |genePredToBed stdin level.s.bed
		awk '$10>1' ../uniq.bed12+ |cut -f 1-12|bedtools intersect -s -wa -a level.s.bed -b - |perl $script_dir/select.pl -i - -r level.gpe -c1 4 -c2 1 -o |grep -i Single |perl -lane '$F[0]=$F[0]."L1";print join "\t",@F;' >>level.gpe.tmp
		awk '$10>1' ../uniq.bed12+ |cut -f 1-12|bedtools intersect -s -wa -a level.s.bed -b - |perl $script_dir/select.pl -i - -r level.gpe -c1 4 -c2 1 |perl -lane '$F[0]=$F[0]."L2";print join "\t",@F;' >>level.gpe.tmp
		export tissue
		perl -lane '$tissue=$ENV{"tissue"};$F[0]=$tissue."-".$F[0];$F[11]=$tissue."-".$F[11];print join "\t",@F;' level.gpe.tmp >level.gpe
		perl $script_dir/grep_exonfasta_from_gpe.pl -r $ref -n -s level.gpe >level.fa
		perl $script_dir/putativeTranslate.pl -i level.fa -f 3 -a
		export cdsLength
		perl -lane 'if ($F[2]-$F[1]+1 <= $ENV{"cdsLength"}) {$F[2]="-";$F[1]="-";}print join "\t",@F'  putativeCDS |perl $script_dir/lncRNA_gpe2putativeCDS_gpe.pl -i - -g level.gpe >level.cds.gpe
		perl $script_dir/grep_cdsfasta_from_gpe.pl -r $ref level.cds.gpe|perl $script_dir/Base.faConvertProtein.fa.pl -i - -o level.cds.fa
		rm final.gtf level.fa long* *.tmp && cd ..
	##2.4 CAGE-seq data processing and supporting information for newly-defined gene models
		###2.4.1 CAGE-seq data processing
		cd /rd1/brick/lisx/pacbio/cage/GR-data/repeat/
		bash run1-cage-GRrawdata.sh
		bash run2_GRctss.sh		
		cd /rd1/brick/lisx/pacbio/cage/FANTOM5/data
		bash run_FANTOM5ctss.sh
		###2.4.2 CAGE-seq supporting informations			
		do cd /rd1/brick/lisx/pacbio/annotation/transcript-cage-${tissue}
 			bash run-cage.sh
		done
				
#3 Gene name assignment (liym@jupiter ~/transcriptome/visualization/geneNameAssign)
	awk '$1!~"Single"' /mnt/share/shenq/forYM/20180702/Final.upgradedMerge2.cds.cageLevel.gpe |genePredToGtf -utr -source=monkeyGene file stdin Final.upgradedMerge2.cds.cageLevel.gtf 
	awk '$1~"Single.level"' /mnt/share/shenq/forYM/20180702/Final.upgradedMerge2.cds.cageLevel.gpe|gpe2bed.pl -t 6 >Final.upgradedMerge2.cds.cageLevel.single.bed6 
	##3.1 Muti-exon genes
		###3.1.1 Assign from macaque refSeq 
		cuffcompare -r /rd1/user/liym/transcriptome/data/rheMac8.refGene.gtf -R -T -o rheMac8.refGene Final.upgradedMerge2.cds.cageLevel.gtf
		awk '$4=="j" || $4=="c" || $4=="="{split($3,a,"|");split($5,b,":");split(b[2],c,"|");print c[2]"\t"a[1]}' rheMac8.refGene.tracking |sort|uniq >inhouse-refGene 
		cuffcompare -r Final.upgradedMerge2.cds.cageLevel.gtf -R -T -o rheMac8.refGene.convert /rd1/user/liym/transcriptome/data/rheMac8.refGene.gtf
		awk '$4=="j" || $4=="c" || $4=="="{split($3,a,"|");split($5,b,":");split(b[2],c,"|");print a[2]"\t"c[1]}' rheMac8.refGene.convert.tracking |sort|uniq >inhouse-refGene.convert
		awk '$9>1' /rd1/user/liym/transcriptome/data/rheMac8.refGene.gpe|cut -f13|join.pl -i1 inhouse-refGene.convert -f1 2 -o1|sort|uniq >tmp;mv tmp inhouse-refGene.convert
		comm -13 inhouse-refGene inhouse-refGene.convert|sort|cut -f1 >refGene.conflict
		join -v 1 -1 1 -2 1 -t $'\t' inhouse-refGene refGene.conflict |awk '{print $0"\tR-refSeq"}' >final_inhouse_geneName.txt
		###3.1.2 Assign from human aligned gencode 
		liftOver -minMatch=0.8 -genePred /rd1/user/liym/transcriptome/data/gencode.v19.annotation.genename.gpe ~/data/liftOver/hg19ToRheMac8.over.chain gencode.v19.TorheMac8.gpe hg19TorheMac8.unmapped
		genePredToGtf -utr file gencode.v19.TorheMac8.gpe gencode.v19.TorheMac8.gtf
		cuffcompare -r gencode.v19.TorheMac8.gtf -R -T Final.upgradedMerge2.cds.cageLevel.gtf
		cuffcompare -r Final.upgradedMerge2.cds.cageLevel.gtf -R -T -o cuffcmp.convert gencode.v19.TorheMac8.gtf
		awk '$4=="=" || $4=="c" || $4=="j"{split($3,a,"|");split($5,b,":");split(b[2],c,"|");print c[2]"\t"a[1]}' cuffcmp.tracking |sort -u >inhouse-humanAlign
		awk '$4=="=" || $4=="c" || $4=="j"{split($3,a,"|");split($5,b,":");split(b[2],c,"|");print a[2]"\t"c[1]}' cuffcmp.convert.tracking |sort -u >inhouse-humanAlign.convert
		awk '$8>1' /rd1/user/liym/transcriptome/data/gencode.v19.annotation.genename.gpe|cut -f12|sort|uniq|join.pl -i1 inhouse-humanAlign.convert -f1 2 -o1|sort|uniq >tmp;mv tmp inhouse-humanAlign.convert
		comm -23 <(cut -f1 final_inhouse_geneName.txt|sort) <(cut -f1 inhouse-humanAlign.convert |sort|uniq -c|awk '$1>1{print $2}'|sort) >inhouse-moreThanOneHgene
		join -v 1 -1 1 -2 1 inhouse-humanAlign inhouse-moreThanOneHgene >inhouse-humanAlign-unique
		comm -23 <(awk '$1!~"Single"' /mnt/share/shenq/forYM/20180702/Final.upgradedMerge2.cds.cageLevel.gpe|cut -f1|sort) <(cut -f1 final_inhouse_geneName.txt|sort)|join.pl -i1 inhouse-humanAlign-unique -o1|awk '{print $1"\t"$2"\tH-Align"}' >>final_inhouse_geneName.txt
		###3.1.3 Assign by those already have gene names
		awk '$1!~"Single"' Final.upgradedMerge2.cds.cageLevel.gpe|awk '{print $12"\t"$1}' >all.geneName.transName
		join.pl -i1 Final.upgradedMerge2.cds.cageLevel.gpe -i2 final_inhouse_geneName.txt -f2 1 -o1|genePredToGtf -utr file stdin assignName.1.gtf
		comm -23 <(cut -f2 all.geneName.transName |sort) <(cat refGene.conflict inhouse-moreThanOneHgene <(cut -f1 final_inhouse_geneName.txt)|sort|uniq)|join.pl -i1 Final.upgradedMerge2.cds.cageLevel.gpe -o1 |genePredToGtf -utr file stdin remainToAssignName.1.gtf
		cuffcompare -r assignName.1.gtf -R -T -o remain1 remainToAssignName.1.gtf
		awk '$4=="c" || $4=="j"' remain1.tracking |awk '{split($3,a,"|");split($5,b,":");split(b[2],c,"|");print a[2]"\t"c[2]}'|join.pl -f1 1 -i2 final_inhouse_geneName.txt -f2 1|cut -f2,4|awk '{print $1"\t"$2"\tAssign-byOther"}' >>final_inhouse_geneName.txt
		###3.1.4 Assign by exon liftover
		comm -23 <(cut -f2 all.geneName.transName|sort) <(cat refGene.conflict inhouse-moreThanOneHgene <(cut -f1 final_inhouse_geneName.txt)|sort)|join.pl -i1 Final.upgradedMerge2.cds.cageLevel.gpe -f1 1 -o1|perl /mnt/share/liym/bin/gpeGetExons.pl >exonForExonLift.rheMac8.bed6
		perl /mnt/share/liym/bin/gpeGetExons.pl -i /rd1/user/liym/transcriptome/data/gencode.v19.annotation.genename.gpe | liftOver -bedPlus=6 -minMatch=0.8 stdin /data/liftover/hg19/hg19ToRheMac8.over.chain.gz gencode.v19.exonTorheMac8.bed6 hg19TorheMac8.exon.unmapped
		grep -v "^#" hg19TorheMac8.unmapped|perl /mnt/share/liym/bin/gpeGetExons.pl|join.pl -i1 gencode.v19.exonTorheMac8.bed6 -f1 4 -f2 4 -o1|awk '{split($4,a,"|");print $1"\t"a[2]}'|sort|uniq|perl /mnt/share/liym/bin/unique.pl -c 2 >exonLift.theSameChr
		awk '{split($4,a,"|");print a[2]"\t"$0}' gencode.v19.exonTorheMac8.bed6|join.pl -i2 exonLift.theSameChr -f2 2 |cut -f2-7 >gencode.v19.exonTorheMac8.noLiftAsWhole.bed6
		bedtools intersect -s -f 1.0 -r -wo -a exonForExonLift.rheMac8.bed6 -b gencode.v19.exonTorheMac8.noLiftAsWhole.bed6 >exonForExonLift.hg19Intersect.bed6+
		awk '{split($4,a,"|");split($10,b,"|");print a[2]"\t"b[1]}' exonForExonLift.hg19Intersect.bed6+ |sort|uniq|perl /mnt/share/liym/bin/unique.pl |awk '{print $1"\t"$2"\tExon-Lift"}' >>final_inhouse_geneName.txt 
		###3.1.5 Assign by those already have gene names again
		comm -23 <(cut -f2 all.geneName.transName |sort) <(cut -f1 final_inhouse_geneName.txt|sort)|sort|uniq|join.pl -i1 Final.upgradedMerge2.cds.cageLevel.gpe -o1 |genePredToGtf -utr file stdin remainToAssignName.2.gtf
		comm -23 <(cut -f1 final_inhouse_geneName.txt|sort) <(awk '{print $12}' assignName.1.gtf|sed 's/;//'|sed 's;";;g'|sort)|join.pl -i1 Final.upgradedMerge2.cds.cageLevel.gpe -o1|genePredToGtf -utr file stdin assignName.2.gtf
		cuffcompare -r assignName.2.gtf -R -T -o remain2 remainToAssignName.2.gtf
		awk '$4=="c" || $4=="j"' remain2.tracking |awk '{split($3,a,"|");split($5,b,":");split(b[2],c,"|");print a[2]"\t"c[2]}'|join.pl -f1 1 -i2 final_inhouse_geneName.txt -f2 1|cut -f2,4|awk '{print $1"\t"$2"\tAssign-byOther2"}' >>final_inhouse_geneName.txt
	#3.2 Single-exon genes
		bedtools intersect -a Final.upgradedMerge2.cds.cageLevel.single.bed6 -b /rd1/user/liym/transcriptome/data/rheMac8.refGene.single.bed6 -wo -s -f 0.2|awk -v OFS="\t" '{split($10,a,"|");print $4,a[1],"R-refSeq"}' >>final_inhouse_geneName.txt
		comm -23 <(cut -f4 Final.upgradedMerge2.cds.cageLevel.single.bed6|sort) <(cut -f1 final_inhouse_geneName.txt|sort)|join.pl -i1 Final.upgradedMerge2.cds.cageLevel.single.bed6 -f1 4 -o1|bedtools intersect -a stdin -b /rd1/user/liym/transcriptome/data/gencode.v19.single.TorheMac8.bed6 -wo -s -f 0.2|awk -v OFS="\t" '{split($10,a,"|");print $4,a[1],"H-Align"}'|sort|uniq|perl /mnt/share/liym/bin/unique.pl -c 1 >>final_inhouse_geneName.txt
		#comm -23 <(cut -f4 Final.upgradedMerge2.cds.cageLevel.single.bed6|sort) <(cut -f1 final_inhouse_geneName.txt|sort)|join.pl -i1 Final.upgradedMerge2.cds.cageLevel.single.bed6 -f1 4 -o1|bedtools intersect -a stdin -b /rd1/user/liym/transcriptome/data/gencode.v19.single.TorheMac8.bed6 -wo -s -f 0.2|awk -v OFS="\t" '{split($10,a,"|");print $4,a[1]}'|sort|uniq|cut -f1|sort|uniq -c|awk '$1>2'|awk 'BEGIN{tag=1}{print $2"\tMutiple-single-"tag"\tMultiple-single";tag+=1;}' >>final_inhouse_geneName.txt 
		#comm -23 <(cut -f4 Final.upgradedMerge2.cds.cageLevel.single.bed6|sort) <(cut -f1 final_inhouse_geneName.txt|sort)|awk 'BEGIN{tag=1}{print $1"\tInhouse-single-"tag"\tInhouse-single";tag+=1}' >>final_inhouse_geneName.txt
	#3.3 Assign by unique locus-geneName pairs
		join.pl -i1 final_inhouse_geneName.txt -i2 transName.geneName.426909|cut -f2,5|sort|uniq|perl /mnt/share/liym/bin/unique.pl -c 2 >XLOC.geneName 
		comm -23 <(awk '$1!~"Single" || ($1~"Single" && $1~"level")' Final.upgradedMerge2.cds.cageLevel.gpe|cut -f1|sort) <(cut -f1 final_inhouse_geneName.txt|sort)|join.pl -i1 transName.geneName.426909 -o1|join.pl -f1 2 -i2 XLOC.geneName -f2 2|awk -v OFS="\t" '{print $1,$3,"LOC_align"}' >>final_inhouse_geneName.txt 
	#3.4 Add inhouse name 
		comm -23 <(cut -f1 transName.geneName.426909|sort) <(cut -f1 final_inhouse_geneName.txt|sort)|join.pl -i1 transName.geneName.426909 -o1|sort -k2,2|awk -v OFS="" -v ORS="" 'BEGIN{tag=1;gene="XLOC_000001"}{if($2==gene){print $1"\t";print "InhouseG";printf("%010d",tag);print "\tInhouse\n";}else{gene=$2;tag+=1;print $1"\t";print "InhouseG";printf("%010d",tag);print "\tInhouse\n";}}' >>final_inhouse_geneName.txt
		#Summary for gene name assignment
		join.pl -i1 final_inhouse_geneName.txt -i2 Final.upgradedMerge2.cds.cageLevel.gpe|cut -f1-3,15|sort -k2,2|awk -v OFS="" -v ORS="" 'BEGIN{gene="7SK";tagG=1;tagT=1;}{if($2!=gene){tagG+=1;gene=$2}print $0"\tMonkeyG";printf("%010d",tagG);print "\tMonkeyT";printf("%010d",tagT);print "\n";tagT+=1}'|awk -v OFS="\t" '{if($3=="Inhouse" || $3=="Multiple"){print $0,$5}else{print $0,$2}}' >geneNameAssign.summary.tsv
		##Header: TCONS inhouse-symbol type XLOC monkeyGeneID monkeyGeneTID finalSymbol
		join.pl -i1 Final.upgradedMerge2.cds.cageLevel.gpe -i2 geneNameAssign.summary.tsv |awk -v OFS="\t" '{print $21,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$22,$13,$14,$15}' >Final.426909.withName.gpe 

#4 The evaluation of newly-defined macaque gene models (@jupiter)
	Ensembl_M=/rd1/user/liym/transcriptome/data/rheMac8.ensGene.96.chr.gpe
	RefSeq_M=/rd1/user/liym/transcriptome/data/rheMac8.refGene.20190613.gpe
	GENCODE_H=/rd1/user/liym/transcriptome/data/gencode.v19.annotation.genename.gpe
	newlyDefined_M=/rd1/user/liym/transcriptome/visualization/Final.426909.withName.gpe
	##3.1 Comparisons with annotations in human dan macaque (Fig1d,e,f)		
		###directory:/rd1/brick/lisx/pacbio/evaluation/Evalu-overlap/
		###with human GENCODE
		liftOver -genePred $newlyDefined_M /share/data/liftover/rheMac8/rheMac8ToHg19.over.chain.gz rheMac8ToHg19.bgm.gpe unmapped.rheMac8ToHg19
		genePredToGtf file rheMac8ToHg19.bgm.gpe rheMac8ToHg19.bgm.gtf
		genePredToGtf file $GENCODE_H hg19.gencode.gtf 
		cuffcompare -r hg19.gencode.gtf -o cuffcompare rheMac8ToHg19.bgm.gtf
		perl /rd1/brick/lisx/scripts/comm_file1_file2.pl -f1 <(less -S hg19.gencode.bed |cut -f 4,13 ) -f2 <(less -S cuffcompare/cuffcompare.tracking |sed 's/gene name/gene_name/g'|awk '{if (($4=="c") || ($4=="j") || ($4=="=") )print $0}'|cut -f 3,5|sed 's/|/\t/g'|sed 's/:/\t/g'|cut -f 2,4,5 ) -n 1 -e1 1 -e2 3|cut -f 2-|sort|uniq > overlap.hg19rheMac8
		###with macaque refSeq
		genePredToGtf file $RefSeq_M rheMac8.refGene.gtf
		genePredToGtf file $newlyDefined_M rheMac8.bgm.gtf 
		cuffcompare -r rheMac8.refGene.gtf -o cuffcompare rheMac8.bgm.gtf
		less -S cuffcompare/cuffcompare.tracking |awk '{if (($4=="c") || ($4=="j") || ($4=="=") )print $0}'|cut -f 3,5|sed 's/|/\t/g'|sed 's/:/\t/g'|cut -f 1,2,4,5|sort|uniq > overlap.refseqBgm
		###with macaque ensembl
		genePredToGtf file $Ensembl_M rheMac8.ensGene.gtf
		cuffcompare -r rheMac8.ensGene.gtf -o cuffcompare_ens rheMac8.bgm.gtf
		less -S cuffcompare/cuffcompare_ens.tracking |awk '{if (($4=="c") || ($4=="j") || ($4=="=") )print $0}'|cut -f 3,5|sed 's/|/\t/g'|sed 's/:/\t/g'|cut -f 1,2,4,5|sort|uniq > overlap.ensBgm
	##3.2 TSS evaluation
		###3.2.1 H3K4me3 data processing (liym@pluto /rd1/user/liym/transcriptome/ChIP-seq; Fig2a)
		perlScript_ym=/rd1/brick/liym/bin
		mv SRR068550.sra brain_pfc.SRR068550.sra
		mv SRR068551.fastq brain_pfc.input.SRR068551.fastq
		mv SRR2927042.fastq brain_ocp.input.SRR2927042.fastq
		mv SRR2927045.fastq brain_put.input.SRR2927045.fastq
		mv SRR2927047.fastq brain_wm.input.SRR2927047.fastq
		mv SRR3030721.fastq brain_cer.input.SRR3030721.fastq
		mv SRR1979002.sra brain_ocp.SRR1979002.sra
		mv SRR1979003.sra brain_wm.SRR1979003.sra
		mv SRR1979004.sra brain_put.SRR1979004.sra
		mv SRR1979005.sra brain_cer.SRR1979005.sra
		$perlScript_ym/fqSeFilter.pl -f 0.1 -b N -q 0.5 -c 20 -a 20 brain_cer.input.SRR3030721.fastq >brain_cer.input.SRR3030721.filter.fastq 2>brain_cer.input.SRR3030721.filter.log 
		$perlScript_ym/fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCG -i brain_cer.SRR1979005.fastq -o brain_cer.SRR1979005.filter.fastq -r brain_cer.SRR1979005.filter.log 
		$perlScript_ym/fqSeFilter.pl -f 0.1 -b N -q 0.5 -c 20 -a 20 brain_ocp.input.SRR2927042.fastq >brain_ocp.input.SRR2927042.filter.fastq 2>brain_ocp.input.SRR2927042.filter.log
		$perlScript_ym/fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCG -i brain_ocp.SRR1979002.fastq -o brain_ocp.SRR1979002.filter.fastq -r brain_ocp.SRR1979002.filter.log 
		$perlScript_ym/fqAdapterFilter.pl -a GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAAAAAA,GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTTAAAAAA -i brain_pfc.input.SRR068551.fastq -o brain_pfc.input.SRR068551.filter.fastq -r brain_pfc.input.SRR068551.filter.log 
		$perlScript_ym/fqSeFilter.pl -f 0.1 -b N -q 0.5 -c 20 -a 20 brain_put.input.SRR2927045.fastq >brain_put.input.SRR2927045.filter.fastq 2>brain_put.input.SRR2927045.filter.log 
		$perlScript_ym/fqAdapterFilter.pl -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCG -i brain_put.SRR1979004.fastq -o brain_put.SRR1979004.filter.fastq -r brain_put.SRR1979004.filter.log 
		$perlScript_ym/fqSeFilter.pl -f 0.1 -b N -q 0.5 -c 20 -a 20 brain_wm.SRR1979003.fastq >brain_wm.SRR1979003.filter.fastq 2>brain_wm.SRR1979003.filter.log 
		ls *.fastq |while read file;do prefix=$(echo $file|sed 's/.fastq//');bwa aln -t 10 /data/bwa/rheMac8/bwtsw $file >${prefix}.out.sai 2>log/${prefix}.aln.log;bwa samse /data/bwa/rheMac8/bwtsw ${prefix}.out.sai $file 2>log/${prefix}.samse.log |samtools view -bSu -| /home/share/local/bin/samtools sort -@ 10 -o ${prefix}.out.sorted.bam - 2>log/${prefix}.samtoolSort.log;done;
		ls *.out.sorted.bam |while read file;do 
			prefix=$(echo $file|sed 's/.out.sorted.bam//');
			bamtools filter -script /mnt/share/liym/data/bwa/bwaAln.filter.json -in $file -out ${prefix}.filtered.sorted.bam;
			bamtools stats -in ${prefix}.filtered.sorted.bam >${prefix}.filtered.sorted.bamStats;
			samtools rmdup -s ${prefix}.filtered.sorted.bam ${prefix}.filtered.rmdup.sorted.bam;
		done;
		ls *rmdup.sorted.bam|grep -v "input"|while read file;do
			prefix=$(echo $file|cut -f1 -d '.');
			factor=$(grep "Total" ${file}Stats|awk '{print 10000000/$3}')
			input=$(grep "Total" ${prefix}.input.*.filter.filtered.rmdup.sorted.bamStats|awk '{print 10000000/$3}')
			tissue=$(echo $prefix|cut -f2 -d '_');
			fragLen=$(grep "predicted fragment length" ${tissue}.macs2.err|awk '{print $13}');
			/mnt/share/share/local/bin/bamCompare -b1 $file -b2 ${prefix}.input.*.filter.filtered.rmdup.sorted.bam --scaleFactors ${factor}:$input --ratio subtract -bs 20 -p 20 --extendReads $fragLen -o ${prefix}.subtract.bw --outFileFormat bigwig >log/${prefix}.bamCompare.log 2>log/${prefix}.bamCompare.err
        done;
		###3.2.2 H3K4me3 signals across TSS (liym@jupiter:/rd1/user/liym/transcriptome/evaluation/TSS_H3K4me3; Fig2b)
		awk -v OFS="\t" '{if($3=="+"){print $2,$4,$4+1,$12,"0",$3}else{print $2,$5-1,$5,$12,"0",$3}}' $newlyDefined_M >Final.426909.withName.TSS.bed6
		awk -v OFS="\t" '{if($3=="+"){print $2,$4,$4+1,$12,"0",$3}else{print $2,$5-1,$5,$12,"0",$3}}' $Ensembl_M >rheMac8.ensGene.96.TSS.bed6
		ls /rd1/user/liym/transcriptome/ChIP-seq/*subtract.bw|while read file;do
			prefix=$(basename $file|sed 's/.subtract.bw//');
			bwtool agg 5000:5000 Final.426909.withName.TSS.bed6,rheMac8.ensGene.96.TSS.bed6 $file ${prefix}.ChIPsignals.inhouseAll.ens.agg.txt;
		done;
		paste *.ChIPsignals.inhouseAll.ens.agg.txt|awk -v OFS="\t" '{print $1,($2+$5+$8+$11+$14)/5,($3+$6+$9+$12+$15)/5}' >ChIPsignals.inhouseAll.ens.5tissueMean.agg.txt
		Rscript ~/bin/lines.R -i=ChIPsignals.inhouseAll.ens.5tissueMean.agg.txt -x="Distance to TSS(bp)" -y="ChIP signals" -c="red,blue" -l="Newly-defined,Ensembl" -y1=0 -y2=4 -o=ChIPsignals.inhouseAll.ens.5tissueMean.agg.pdf 
		###3.2.3 CpG density across TSS
		wget http://hgdownload.soe.ucsc.edu/goldenPath/rheMac8/database/cpgIslandExt.txt.gz .
		cut -f 2-5 cpgIslandExt.txt >CpG.bed4
		cat $newlyDefined_M|perl -lane '$a[0]=$F[1];if ($F[2] eq "+"){$a[1]=$F[3];$a[2]=$F[3]+1;} else {$a[1]=$F[4]-1;$a[2]=$F[4];} $a[3]=$F[0];$a[4]=$F[7];$a[5]=$F[2];print join "\t",@a;'|sort -k1,1 -k2,2n|bedtools closest -D "a" -t "first" -a - -b <(sort -k1,1 -k2,2n CpG.bed4) > dis2cpg.all
		cat $Ensembl_M |cut -f 2-|perl -lane '$a[0]=$F[1];if ($F[2] eq "+"){$a[1]=$F[3];$a[2]=$F[3]+1;} else {$a[1]=$F[4]-1;$a[2]=$F[4];} $a[3]=$F[0];$a[4]=$F[7];$a[5]=$F[2];print join "\t",@a;'|sort -k1,1 -k2,2n|bedtools closest -D "a" -t "first" -a - -b <(sort -k1,1 -k2,2n CpG.bed4) > dis2cpg.ensgene
	##3.3 TTS evaluation (liym@jupiter ~/transcriptome/evaluation/polyASeq; Fig2c)
		#download UTL:wget http://hgdownload.soe.ucsc.edu/gbdb/rheMac2/bbi/polyASeqSites*Fwd.bw 
					 #wget http://hgdownload.soe.ucsc.edu/gbdb/rheMac2/bbi/polyASeqSites*Rev.bw
		cat rawBW/*Fwd.bw.bg |perl -lane '$n++;$F[4]=$F[3];$F[3]="PolyA_".$n;$F[5]="+";print join "\t",@F;' > polyA.rheMac2.bed6 
		cat rawBW/*Rev.bw.bg | perl -lane '$n++;$F[4]=$F[3];$F[3]="PolyA_".$n;$F[5]="-";print join "\t",@F;' >> polyA.rheMac2.bed6
		liftOver polyA.rheMac2.bed6 /mnt/share/share/data/liftover/rheMac2/rheMac2ToRheMac8.over.chain.gz polyA.rheMac8.bed6 unmap 
		cat $newlyDefined_M|perl -lane '$a[0]=$F[1];if ($F[2] eq "+"){$a[1]=$F[4]-1;$a[2]=$F[4];} else {$a[1]=$F[3];$a[2]=$F[3]+1;} $a[3]=$F[0];$a[4]=$F[7];$a[5]=$F[2];print join "\t",@a;'|sort -k1,1 -k2,2n|bedtools closest -s -D "a" -t "first" -a - -b <(sort -k1,1 -k2,2n GSE30198.rheMac8.polyASeq.bed6)  >dis2polya.all
		cat $Ensembl_M|perl -lane '$a[0]=$F[1];if ($F[2] eq "+"){$a[1]=$F[4]-1;$a[2]=$F[4];} else {$a[1]=$F[3];$a[2]=$F[3]+1;} $a[3]=$F[0];$a[4]=$F[7];$a[5]=$F[2];print join "\t",@a;'|sort -k1,1 -k2,2n|bedtools closest -s -D "a" -t "first" -a - -b <(sort -k1,1 -k2,2n GSE30198.rheMac8.polyASeq.bed6) >dis2polya.ensgene
		cat <(awk -v OFS="\t" '{print $13,"Newly-defined"}' dis2polya.all) <(awk -v OFS="\t" '{print $13,"Ensembl"}' dis2polya.ensgene) >dis2ploya.inhouse.ensembl
	##3.4 Splice site evaluation (liym@jupiter ~/transcriptome/evaluation/ss_motif;Fig2d)
		perl ~/bin/gpeFeature.pl -i ~/transcriptome/visualization/Final.426909.withName.gpe|awk -v OFS="\t" '{if($6=="+"){print $1,$2-3,$2+6,$4,$5,$6}else{print $1,$3-6,$3+3,$4,$5,$6}}' | fastaFromBed -name -s -fi ~/data/genome/rheMac8/rheMac8.fa -bed stdin -fo stdout|grep -v "^>" - >Final.426909.5SS.motif
		perl ~/bin/gpeFeature.pl -i ~/transcriptome/visualization/Final.426909.withName.gpe|awk -v OFS="\t" '{if($6=="+"){print $1,$3-20,$3+3,$4,$5,$6}else{print $1,$2-3,$2+20,$4,$5,$6}}' | fastaFromBed -name -s -fi ~/data/genome/rheMac8/rheMac8.fa -bed stdin -fo stdout|grep -v "^>" - >Final.426909.3SS.motif
		python /mnt/share/liym/tools/weblogo-master/weblogo -f Final.426909.5SS.motif -o 5ss.weblogo.pdf -F pdf -A dna -t "5ss" --number-interval 1 --fineprint "" --errorbars False -C green A 'Adenine' -C orange G 'guanine' -C red T 'Thymine' -C blue C 'Cytosine'
		bin/python /mnt/share/liym/tools/weblogo-master/weblogo -f Final.426909.3SS.motif -o 3ss.weblogo.pdf -F pdf -A dna -t "3ss" --number-interval 1 --fineprint "" --errorbars False -C green A 'Adenine' -C orange G 'guanine' -C red T 'Thymine' -C blue C 'Cytosine'
	##3.5 Protein sequence covergae (lixs@pluto:/rd1/brick/lisx/pacbio/NCBI_inputdata/protein-ORF/; Fig2e)
		perlScript_lsx=/rd1/brick/lisx/scripts
		perl $perlScript_lsx/bin-mergeValue.pl -f <(perl /rd1/brick/lisx/scripts/gpe-utr-cds-length.pl -f1 hg19.refGene.20190613.gpe |cut -f 1,14|sed 's/:/\t/g'|awk '{if ($2!=$5)print $5"\t"$2 }' ) |sed 's/,/\t/g'|awk '{a=$2;for (i=2;i<=NF;i++){if (a<$i){a=$i}};print $1"\t"a  }' > hg19.cdslength
		perl $perlScript_lsx/bin-mergeValue.pl -f <(perl /rd1/brick/lisx/scripts/gpe-utr-cds-length.pl -f1 <(awk '{if ($6!=$7){print "0\t"$0}}' rheMac8.bgm.gpe )|cut -f 1,14|sed 's/:/\t/g'|awk '{print $5"\t"$2 }')|sed 's/,/\t/g'|awk '{a=$2;for (i=2;i<=NF;i++){if (a<$i){a=$i}};print $1"\t"a  }'  > rheMac8.cdslength
		perl $perlScript_lsx/comm_file1_file2.pl -f1 hg19.cdslength -f2 rheMac8.cdslength -n 1 -e1 1 -e2 1 |cut -f 1,2,4 > hg19RheMac8.cdslength
		### run blastp for bgm transcriptome 
		cd blastp
		bash run-blastp.sh 
		### process blastp result and statistic the result
		cd ../blastp-orf
		bash run-blRes-CovIden.sh all 75 
		bash run-blCovIden-statistic.sh all.iden75.Topres
		bash run-finalfigure75.sh
	##3.6 CDS comparisions with human CDS(lisx@pluto:/rd1/brick/lisx/pacbio/evaluation/cds/)
		###3.6.1 CDS length comparison with human refSeq annotations (Fig2f)
			perl $perlScript_lsx/bin-mergeValue.pl -f <(perl ${perlPath}/gpe-utr-cds-length.pl -f1 <(awk '{print $0}' /share/data/structure/gpe/hg19.refGene.gpe ) |cut -f 1,14|sed 's/:/\t/g'|awk '{if ($2!=$5)print $5"\t"$2 }' ) |sed 's/,/\t/g'|awk '{a=$2;for (i=2;i<=NF;i++){if (a<$i){a=$i}};print $1"\t"a  }' > hg19.cdslength
			perl $perlScript_lsx/bin-mergeValue.pl -f <(perl ${perlPath}/gpe-utr-cds-length.pl -f1 <(awk '{if ($6!=$7){print "0\t"$0}}' rheMac8.bgm.gpe )|cut -f 1,14|sed 's/:/\t/g'|awk '{print $5"\t"$2 }')|sed 's/,/\t/g'|awk '{a=$2;for (i=2;i<=NF;i++){if (a<$i){a=$i}};print $1"\t"a  }'  > rheMac8.cdslength
			perl $perlScript_lsx/comm_file1_file2.pl -f1 hg19.cdslength -f2 rheMac8.cdslength -n 1 -e1 1 -e2 1 |cut -f 1,2,4 > hg19RheMac8.cdslength
			perl $perlScript_lsx/bin_value.pl -f <(less -S hg19RheMac8.cdslength |awk '{print $1"\t"($3-$2)/$2+0.05}' ) -value 2 -n 40 -max 3 -min -1 -equal max -outbin min~max -outnumber percent|awk '{a=a+$2;print $0"\t"a}'> rhe-hg.box
		###3.6.2 Coverage proportation and identity (Fig2g,h)
			cd AAcompare/
			### get the homologous genes of human and macaque
			perl $perlScript_lsx/comm_file1_file2.pl -f1 <(less -S ../hg19RheMac8.cdslength |cut -f 1 ) -f2 <(less -S /share/data/structure/gpe/hg19.refGene.gpe |awk '{print $13"\t"$0}') -n 1 -e1 0 -e2 20|cut -f 4- > hg19.refHomo.gpe
			perl $perlScript_lsx/comm_file1_file2.pl -f1 <(less -S ../hg19RheMac8.cdslength |cut -f 1 ) -f2 <(less -S ../rheMac8.bgm.gpe |awk '{print $12"\t"$0}') -n 1 -e1 0 -e2 20|cut -f 3- > rheMac8.bgmHomo.gpe
			### get the aa sequence of human and macaque homologous genes
			#### for macaque
			cd rheMac8-aa
			less ../rheMac8.bgmHomo.gpe|grep -v chrUn|awk '{if ($6!=$7)print $0}' > rheMac8.gpe
			perl ${perlPath}/gpe-outputCDS.pl -f1 <(awk '{print 0"\t"$0}' rheMac8.gpe) > rheMac8.gpe.cds
			gpe2bed.pl <(cut -f 2- rheMac8.gpe.cds) > rheMac8.bed
			bedtools getfasta -fi /share/data/fna/rheMac8/all.fa -bed rheMac8.bed -fo rheMac8.fa -split -tab -s
			perl ${perlPath}/comm_file1_file2.pl -f1 <(less -S rheMac8.gpe|awk '{print ">"$1"\t>"$1"#"$12}') -f2 <(paste <(cut -f 1-4 rheMac8.bed) rheMac8.fa |awk '{print ">"$4"\t"$6}' ) -n 1 -e1 1 -e2 1|cut -f 2,4|sed 's/\t/\n/g' > trans.fa
			perl ${perlPath}/Base.faConvertProtein.fa.pl --Fasta trans.fa --output trans.aa
			makeblastdb -in trans.aa -dbtype prot -parse_seqids -out trans.aa > blastdb.log
			rm rheMac8.*
			less -S trans.aa|awk '{if(0==NR%2)printf("%s\n",$0);else printf("%s\t",$0)}'|awk '{print $1"\t"length($2)}'|sed 's/>//g' > trans.aa.length
			#### for human
			cd ../hg19-aa
			less -S ../hg19.refHomo.gpe |awk '{if ($2!~/_/)print $0}'|awk '{if ($6!=$7)print $0}' > hg19.gpe
			perl ${perlPath}/gpe-outputCDS.pl -f1 <(awk '{print 0"\t"$0}' hg19.gpe) > hg19.gpe.cds
			gpe2bed.pl <(cut -f 2- hg19.gpe.cds) > hg19.bed
			bedtools getfasta -fi /share/data/fna/hg19/all.fa -bed hg19.bed -fo hg19.fa -split -tab -s
			perl ${perlPath}/comm_file1_file2.pl -f1 <(less -S hg19.gpe|awk '{print ">"$1"\t>"$1"#"$12}') -f2 <(paste <(cut -f 1-4 hg19.bed) hg19.fa |awk '{print ">"$4"\t"$6}' ) -n 1 -e1 1 -e2 1|cut -f 2,4|sed 's/\t/\n/g' > trans.fa
			perl ${perlPath}/Base.faConvertProtein.fa.pl --Fasta trans.fa --output trans.aa 2>conPro.log
			rm hg19.*
			less -S trans.aa|awk '{if(0==NR%2)printf("%s\n",$0);else printf("%s\t",$0)}'|awk '{print $1"\t"length($2)}'|sed 's/>//g' > trans.aa.length
			### blastp
			cd ../result
			blastp -query hg19-aa/trans.aa -out result/blastp.out -db rheMac8-aa/trans.aa -outfmt 6 -evalue 1e-5 -num_threads 50 > result/blastp.err 2> result/blastp.log
			### calculate the identity and coverage
			cd ../similarity
			less -S ../result/blastp.out |sed 's/#/\t/g'|awk '{if ($2==$4)print $0}' >blastp.out
			perl $perlScript_lsx/comm_file1_file2.pl -f1 ../rheMac8-aa/trans.aa.length -f2 <(perl $perlScript_lsx/comm_file1_file2.pl -f1 ../hg19-aa/trans.aa.length -f2 blastp.out -n 1 -e1 1 -e2 1 |awk '{print $4"\t"$0}') -n 1 -e1 1 -e2 1 |cut -f 1,2,4,5,8-|sed 's/#/\t/g'|awk '{print $2"\t"$0}'|awk '{if ($1==$6)print $0}'|cut -f 1,2,4,5,7- > blastp.result
			perl $perlScript_lsx/bin_value.pl -f <(cut -f 6 blastp.result) -value 1 -n 20 -max 100 -min 0 -equal max -outbin min~max -outnumber percent |awk '{a=a+$2;print $0"\t"a}'>distri.allidentity
			less -S blastp.result |awk '{print $1"\t"($11-$10+1)/$5"\t"$6}'|sort|uniq > all.coverage.identity
			perl $perlScript_lsx/bin_value.pl -f <(awk '{print $2*100}' all.coverage.identity ) -value 1 -n 20 -max 100 -min 0 -equal max -outbin min~max -outnumber percent |awk '{a=a+$2;print $0"\t"a}' >distri.allcoverage
			perl $perlScript_lsx/bin-mergeValue.pl -f <(cut -f 1,3 all.coverage.identity|sort|uniq ) |sed 's/,/\t/g'|awk '{a=$2;for (i=2;i<=NF;i++){if (a<$i){a=$i}};print $1"\t"a  }' > all.topidnetity
			perl $perlScript_lsx/bin_value.pl -f all.topidnetity -value 2 -n 20 -max 100 -min 0 -equal max -outbin min~max -outnumber percent |awk '{a=a+$2;print $0"\t"a}' > distri.topidentity
			perl $perlScript_lsx/bin-mergeValue.pl -f <(perl $perlScript_lsx/comm_file1_file2.pl -f1 all.topidnetity -f2 <(awk '{print $1"\t"$3"\t"$2}' all.coverage.identity ) -n 2 -e1 0 -e2 1|awk '{print $3"\t"$5*100}' ) |sed 's/,/\t/g'|awk '{a=$2;for (i=2;i<=NF;i++){if (a<$i){a=$i}};print $1"\t"a  }' > all.topidnetity-topcoverage
			perl $perlScript_lsx/bin_value.pl -f all.topidnetity-topcoverage  -value 2 -n 20 -max 100 -min 0 -equal max -outbin min~max -outnumber percent |awk '{a=a+$2;print $0"\t"a}' > distri.topidnetity-topcoverage

#5 PA identification and evaluation(liym@jupiter ~/transcriptome/PA)
	##5.1 PA usage on genes in human and macaque tissues 
		ls ../hg19PacBio/*bed12+|grep -v "hCere"|while read file;do
			tissue=$(basename $file|sed 's/.processed.sky.bed12+//');
			sh ~/transcriptome/scripts/PAratio.OnGenes.sh -g ~/transcriptome/data/gencode.v19.annotation.genename.gpe -r $file -o $tissue &
		done;
		ls /mnt/share/liym/transcriptome/pacBioData/rheMac8/*pass1.bed12+|while read file;do
			tissue=$(basename $file|cut -f1 -d '.');
			sh ~/transcriptome/scripts/PAratio.OnGenes.sh -g ~/transcriptome/visualization/Final.426909.withName.gpe -r $file -o $tissue &
		done;
		ls */PA.onGene.tsv |while read file;do
			prefix=$(dirname $file);  
        	 -v OFS="\t" '{split($4,a,",");split($5,readN,",");split($6,usage,",");for(i=1;i<=length(a);i++){print $1,a[i]-1,a[i],$3,readN[i],$2,usage[i]}}' $file >${prefix}.PAusage.bed6+;
        done;
	##5.2 PA evaluation in macaque tissues
		mkdir evaluation & cd evaluation
		ln -s ../brain.PAusage.bed6+ rheMac8.brain.PAusage.bed6+
		ln -s ../cerebellum.PAusage.bed6+ rheMac8.cerebellum.PAusage.bed6+
		ln -s ../heart.PAusage.bed6+ rheMac8.heart.PAusage.bed6+
		ln -s ../testis.PAusage.bed6+ rheMac8.testis.PAusage.bed6+
		###5.2.1 Nucleotides distribution across point (Figure S3)
			ls rheMac8*PAusage.bed6+|while read file;do
				prefix=$(echo $file|cut -f1,2 -d '.');
				awk '$5>1' $file|perl /mnt/share/liym/bin/nucleotide_count_AcrossPoint.pl -u -100 -d 100 -f /mnt/share/liym/data/genome/rheMac8/rheMac8.fa -s >${prefix}.Gt1.nucleotide.tsv 2>>nucleostideDistribution.Gt1.log
				grep -v "Position" ${prefix}.Gt1.nucleotide.tsv |awk -v OFS="\t" '{print $1,$2/($2+$3+$4+$5),$3/($2+$3+$4+$5),$4/($2+$3+$4+$5),$5/($2+$3+$4+$5)}'|sort -k1 -n|Rscript /mnt/share/liym/bin/lines.R -y1=0 -y2=0.6 -x="Distance to cleavage site(bp)" -y="Frequency" -c="red,blue,green,purple" -l="A,C,G,U" -o=${prefix}.Gt1.nucleotide.pdf
			done;
		###5.2.2 PA motif distribution (Figure S2)
			ls rheMac8*PAusage.bed6+|while read file;do
				prefix=$(echo $file|cut -f1,2 -d '.');
				awk '$5>1' $file|PAmotifFreq.agg.pl -m /mnt/share/liym/data/PAmotif.12.txt -u 50 -f /mnt/share/liym/data/genome/rheMac8/rheMac8.fa >${prefix}.PAmotif.Gt1.agg.txt;
				grep -v "pos" ${prefix}.PAmotif.Gt1.agg.txt|Rscript /mnt/share/liym/bin/lines.R -l=AAUAAA,AAUACA,AAUAGA,AAUAUA,AAUGAA,ACUAAA,AGUAAA,AUUAAA,CAUAAA,GAUAAA,UAUAAA -o=${prefix}.PAmotif.Gt1.agg.pdf -y1=0 -y2=0.06 -x="Distance to cleavage site(bp)" -y="Frequency" -c=#FF0000FF,#FF8000FF,#FFFF00FF,#80FF00FF,#00FF00FF,#00FF80FF,#00FFFFFF,#0080FFFF,#0000FFFF,#8000FFFF,#FF00FFFF,#FF0080FF 
			done;
		###5.2.3 PA usage vs. PA motif (Fig.4a & Figure S5)
			ls rheMac8*bed6+|sed 's/.bed6+//'|while read file;do
				perl /mnt/share/liym/bin/PAmotifFinding.pl -b ${file}.bed6+ -l 50 -f /mnt/share/liym/data/genome/rheMac8/rheMac8.fa >${file}.PAmotif.bed6+ &
			done;
		###5.2.4 Compare with poly(A)-seq (Figure S4)
			cd ../PA-seq
			awk -v OFS="\t" -v total=2615605 '{count=sprintf("%0.f",$5*total/1000000);print $1,$2,$3,$4,count,$6}' polyASeqSitesBrain.RPM.rheMac8.bed6 >polyASeqSitesBrain.count.bed6 
			awk -v OFS="\t" -v total=4836387 '{count=sprintf("%0.f",$5*total/1000000);print $1,$2,$3,$4,count,$6}' polyASeqSitesTestis.RPM.rheMac8.bed6 >polyASeqSitesTestis.count.bed6 
			awk -v OFS="\t" '{print $2,$4,$5,$12,"0",$3}' ~/transcriptome/visualization/Final.426909.withName.gpe|sort|uniq|bedtools intersect -a polyASeqSitesBrain.count.bed6 -b stdin -wo |cut -f1-6,10 |sort|uniq >brain.PA.intersectGene.tsv
			awk -v OFS="\t" '{print $2,$4,$5,$12,"0",$3}' ~/transcriptome/visualization/Final.426909.withName.gpe|sort|uniq|bedtools intersect -a polyASeqSitesTestis.count.bed6 -b stdin -wo |cut -f1-6,10 |sort|uniq >testis.PA.intersectGene.tsv
			tissue=brain
			firstGname=$(perl /mnt/share/liym/bin/unique.pl -i ${tissue}.PA.intersectGene.tsv -c 4|sort -k7,7|head -n1|cut -f7)
			perl /mnt/share/liym/bin/unique.pl -i ${tissue}.PA.intersectGene.tsv -c 4|sort -k7,7|awk -v OFS="\t" -v gene=$firstGname '
					BEGIN{geneName=gene;PA="";readNum="";sum=0;chr="";strand=""}
					{
							if($7 == geneName){
									if(sum==0){
											geneName=$7;chr=$1;strand=$6;sum+=$5;PA=$3;readNum=$5;
									}else{
											sum+=$5;PA=PA","$3;readNum=readNum","$5;
									}
							}else{
									if(readNum~","){
											split(readNum,a,",");ratio=sprintf("%.4f",a[1]/sum);for(i=2;i<=length(a);i++){ratio=ratio","sprintf("%.4f",a[i]/sum)}
									}else{
											ratio=1.00 
									}
							print chr,strand,geneName,PA,readNum,ratio;
							geneName=$7;chr=$1;strand=$6;sum=$5;PA=$3;readNum=$5;
							}
					}
					END{
							if(readNum~","){
									split(readNum,a,",");ratio=sprintf("%.4f",a[1]/sum);for(i=2;i<=length(a);i++){ratio=ratio","sprintf("%.4f",a[i]/sum)}
							}else{
									ratio=1.00 
							}
							print chr,strand,geneName,PA,readNum,ratio;
			
					}' >${tissue}.PA.onGene.tsv 
			ls *PA.onGene.tsv |while read file;do
				prefix=$(echo $file|cut -f1 -d '.');  
				awk -v OFS="\t" '{split($4,a,",");split($5,readN,",");split($6,usage,",");for(i=1;i<=length(a);i++){print $1,a[i]-1,a[i],$3,readN[i],$2,usage[i]}}' $file >${prefix}.PAusage.bed6+;
			done;   
			cd ../evaluation;
			awk -v OFS="\t" '{if($2<30){print $1,"0",$3+30,$4,$5,$6,$7}else{print $1,$2-30,$3+30,$4,$5,$6,$7}}' rheMac8.brain.PAusage.bed6+|bedtools intersect -a stdin -b ../PA-seq/rheMac8/brain.PAusage.bed6+ -wo >brain.PA.isoSeq.PAseq.usage.tsv
			awk -v OFS="\t" '{if($2<30){print $1,"0",$3+30,$4,$5,$6,$7}else{print $1,$2-30,$3+30,$4,$5,$6,$7}}' rheMac8.testis.PAusage.bed6+|bedtools intersect -a stdin -b ../PA-seq/rheMac8/testis.PAusage.bed6+ -wo >testis.PA.isoSeq.PAseq.usage.tsv

#6 Cis- and trans-regulations of PA sites selection
	##6.1 Distal and proximal PA strength comparison (Fig.4b & Figure S6)
		ls ../*/PA.onGene.tsv|grep -vE "h[A-Z]|mBrain" |while read file;do
			prefix=$(echo $file|cut -f2 -d '/');  
			awk -v OFS="\t" '$4~","{split($4,a,",");split($5,readN,",");split($6,usage,",");if($2=="+"){print $1,a[1]-1,a[1],$3,readN[1],$2,usage[1]}else{print $1,a[length(a)]-1,a[length(a)],$3,readN[length(a)],$2,usage[length(a)]}}' $file >${prefix}.5most.PAusage.bed6+;
			awk -v OFS="\t" '$4~","{split($4,a,",");split($5,readN,",");split($6,usage,",");if($2=="-"){print $1,a[1]-1,a[1],$3,readN[1],$2,usage[1]}else{print $1,a[length(a)]-1,a[length(a)],$3,readN[length(a)],$2,usage[length(a)]}}' $file >${prefix}.3most.PAusage.bed6+;
		done;
		ls *most*bed6+|while read file;do prefix=$(echo $file|cut -f1,2 -d '.');perl /mnt/share/liym/bin/PAmotifFinding.pl -b $file -l 50 -f /mnt/share/liym/data/genome/rheMac8/rheMac8.fa >${prefix}.PAmotif.tsv;done;
		ls *3most.PAmotif.tsv|while read file;do
			tissue=$(echo $file|cut -f1 -d '.');
			paste ${tissue}.5most.PAmotif.tsv ${tissue}.3most.PAmotif.tsv |awk '$5>1 && $16>1'|cut -f7,8,18,19|awk -v OFS="\t" '{if($2=="AATAAA"){print "1-"$2,$3}else if($2=="NA"){print "4-No",$3}else if($2=="ATTAAA"){print "2-"$2,$3}else{print "3-Other",$3}}' |Rscript /mnt/share/liym/bin/boxplot_formula.R -y1=0 -y2=1 -x="Motif of proximal PA" -y="PA usage" -outline=F -pvalues=T -o=${tissue}.3mostPAusage.binnedByPAS.boxplot.pdf
			paste ${tissue}.5most.PAmotif.tsv ${tissue}.3most.PAmotif.tsv |awk '$5>1 && $16>1'|cut -f7,8,18,19|awk -v OFS="\t" '{if($2=="AATAAA"){print "1-"$2,$3}else if($2=="NA"){print "3-No",$3}else{print "2-Other",$3}}' |Rscript /mnt/share/liym/bin/boxplot_formula.R -y1=0 -y2=1 -x="Motif of proximal PA" -y="PA usage" -outline=F -pvalue=T -o=${tissue}.3mostPAusage.binnedByPAS.3bin.boxplot.pdf
        done;
	##6.2 3D hist for distal PA sites (Fig.4c & Figure S7)
		ls *3most.PAmotif.tsv|while read file;do
			tissue=$(echo $file|cut -f1 -d '.');
			paste ${tissue}.5most.PAmotif.tsv ${tissue}.3most.PAmotif.tsv |awk '$5>1 && $16>1'|cut -f7,8,18,19|awk -v OFS="\t" '{if($2!="AATAAA" && $2!="NA"){if($4!="AATAAA" && $4!="NA"){print "5-Other","3-Other",$3}else{print "5-Other","3-"$4,$3}}else{if($4!="AATAAA" && $4!="NA"){print "5-"$2,"3-Other",$3}else{print "5-"$2,"3-"$4,$3}}}' >${tissue}.3mostPAusage.binnedByPAS.3bin.tsv
        done;
	##6.3 3'UTR length comparison between brain and testis (lixs@pluto:/rd1/brick/lisx/pacbio/APA; Fig.4d)
		### calculate the weighted 3'UTR length for tissues (case:brain, the same with heart,liver,testis )
		cd brain
		bash run-3UTRlength.sh
		cd ../
		### pairwise comparision of 3'UTR length
		bash run-newcompare.sh
	##6.4 PA usage comparison between species and tissues (liym@jupiter:~/transcriptome/PA; Fig.4e)
		ls h*PAusage.bed6+|grep -v "To"|grep -v "heart.PAusage.bed6+"|sed 's/.PAusage.bed6+//'|while read file;do liftOver -bedPlus=6 ${file}.PAusage.bed6+ /data/liftover/hg19/hg19ToRheMac8.over.chain.gz ${file}.TorheMac8.PAusage.bed6+ ${file}.unmapped & done;
		bedtools closest -s -d -a <(awk '$5>1' brain.PAusage.bed6+) -b <(awk '$5>1' cerebellum.PAusage.bed6+)|bedtools closest -s -d -a stdin -b <(awk '$5>1' heart.PAusage.bed6+) |bedtools closest -s -d -a stdin -b <(awk '$5>1' testis.PAusage.bed6+)|awk '$15<=30 || $23<=30 || $31<=30' >Rbcht.PAclosest.usage.PAGt1.tsv
		bedtools closest -s -d -a <(awk '$5>1' hBrain.TorheMac8.PAusage.bed6+) -b <(awk '$5>1' hCere.TorheMac8.PAusage.bed6+)|bedtools closest -s -d -a stdin -b <(awk '$5>1' hHeart.TorheMac8.PAusage.bed6+) |bedtools closest -s -d -a stdin -b <(awk '$5>1' hLiver.TorheMac8.PAusage.bed6+)|awk '$15<=30 || $23<=30 || $31<=30' >Hbchl.PAclosest.usage.PAGt1.tsv
		ls *PAusage.bed6+|grep -vE "To|old|^m|^h[A-Z]"|while read file;do
			tissue=$(echo $file|cut -f1 -d '.');
			bedtools closest -s -d -a <(awk '$5>1' $file) -b <(awk '$5>1' hBrain.TorheMac8.PAusage.bed6+) |bedtools closest -s -d -a stdin -b <(awk '$5>1' hCere.TorheMac8.PAusage.bed6+) |bedtools closest -s -d -a stdin -b <(awk '$5>1' hHeart.TorheMac8.PAusage.bed6+) |bedtools closest -s -d -a stdin -b <(awk '$5>1' hLiver.TorheMac8.PAusage.bed6+) |awk '$15<=30 || $23<=30 || $31<=30 || $39<=30'|cut -f1-6 >R.${tissue}.PAclosest.usage.Gt1.bed6
        done;
		cat *.PAclosest.usage.Gt1.bed6|awk -v OFS="\t" '{print $1,$2,$3,$4,"0",$6}'|sort -k6,6 -k1,1 -k2,2n|uniq|awk -v OFS="\t" '
				BEGIN{chr="chr1";start=63567;end=63568;name="KLHL17";strand="+";}
				{
						if($1==chr && $4==name && $6==strand){
								if($3-end<=30){end=$3}else{print chr,start,end,name,"0",strand;start=$2;end=$3}
						}else{
								print chr,start,end,name,"0",strand;
								chr=$1;start=$2;end=$3;name=$4;strand=$6;
						}
				}
				END{print chr,start,end,name,"0",strand;}
		'| bedtools closest -s -d -a stdin -b <(awk '$5>1' brain.PAusage.bed6+)|bedtools closest -s -d -a stdin -b <(awk '$5>1' cerebellum.PAusage.bed6+)|bedtools closest -s -d -a stdin -b <(awk '$5>1' heart.PAusage.bed6+) |bedtools closest -s -d -a stdin -b <(awk '$5>1' testis.PAusage.bed6+)| bedtools closest -s -d -a stdin -b <(awk '$5>1' hBrain.TorheMac8.PAusage.bed6+) |bedtools closest -s -d -a stdin -b <(awk '$5>1' hCere.TorheMac8.PAusage.bed6+) |bedtools closest -s -d -a stdin -b <(awk '$5>1' hHeart.TorheMac8.PAusage.bed6+) |bedtools closest -s -d -a stdin -b <(awk '$5>1' hLiver.TorheMac8.PAusage.bed6+)  >Rbcht.Hbchl.PAclosest.usage.PAGt1.tsv
		awk -v OFS="" -v ORS="" -v FS="\t" '{
				print $1":"$6":"$2"-"$3;
				for(i=14;i<=NF;i+=8){
						if($i<=30){
								print "\t"$(i-1);
						}else{
								print "\tNA";
						}
				}
				print "\n";
		}' Rbcht.Hbchl.PAclosest.usage.PAGt1.tsv|sort|uniq >Rbcht.Hbchl.PAclosest.usage.PAGt1.final.tsv

#7 Identification of human-specific distal PA sites and transcriptome evolution(liym@jupiter ~/transcriptome/PA/crossSpecies/speciesSpecific/version_all)
	##7.1 Identification of human-specific distal PA sites
		###7.1.1 Only keep gene with all PA sites lifted over to rheMac8
			awk -v OFS="\t" '{if($4!~","){print $1,$4-1,$4,$3":"$4,$5,$2}else{split($4,a,",");split($5,b,",");for(i=1;i<=length(a);i++){print $1,a[i]-1,a[i],$3":"a[i],b[i],$2}}}' ../../../hBrain/PA.onGene.tsv >H.brain.PA.onGene.bed6 
			awk -v OFS="\t" '{if($4!~","){print $1,$4-1,$4,$3":"$4,$5,$2}else{split($4,a,",");split($5,b,",");for(i=1;i<=length(a);i++){print $1,a[i]-1,a[i],$3":"a[i],b[i],$2}}}' ../../../hCere/PA.onGene.tsv >H.cerebellum.PA.onGene.bed6 
			awk -v OFS="\t" '{if($4!~","){print $1,$4-1,$4,$3":"$4,$5,$2}else{split($4,a,",");split($5,b,",");for(i=1;i<=length(a);i++){print $1,a[i]-1,a[i],$3":"a[i],b[i],$2}}}' ../../../hHeart/PA.onGene.tsv >H.heart.PA.onGene.bed6 
			ls H.*PA.onGene.bed6|while read file;do
				tissue=$(echo $file|cut -f2 -d '.');
				liftOver -minMatch=0.75 H.${tissue}.PA.onGene.bed6 /data/liftover/hg19/hg19ToRheMac8.over.chain.gz H.${tissue}.PA.onGene.ToRheMac8.bed6 H.${tissue}.ToRheMac8.unmapped &
			done;
			ls *PA.onGene.ToRheMac8.bed6|sed 's/.bed6//'|while read file;do
				tissue=$(echo $file|cut -f2 -d '.');
				liftOver -minMatch=0.75 ${file}.bed6 /data/liftover/rheMac8/rheMac8ToHg19.over.chain.gz ${file}.Tohg19.bed6 ${tissue}.ToRheMac8.Tohg19.unmapped &
			done;
			ls H.*.PA.onGene.ToRheMac8.bed6|while read file;do
				tissue=$(basename $file|cut -f2 -d '.');
				join.pl -i1 H.${tissue}.PA.onGene.bed6 -f1 4 -i2 H.${tissue}.PA.onGene.ToRheMac8.Tohg19.bed6 -f2 4 |awk '$1==$7 && $2==$8 && $6==$12'|cut -f4|join.pl -i1 H.${tissue}.PA.onGene.ToRheMac8.bed6 -f1 4 -o1|awk -v OFS="\t" '{split($4,a,":");print a[1],$0}' |join.pl -i2 <(cat <(grep -v "#" *${tissue}.ToRheMac8.*unmapped|cut -f4 |cut -f1 -d ":") <(join.pl -i1 H.${tissue}.PA.onGene.bed6 -f1 4 -i2 H.${tissue}.PA.onGene.ToRheMac8.Tohg19.bed6 -f2 4 |awk '$1!=$7 || $2!=$8 || $6!=$12'|cut -f4 |cut -f1 -d ":")|sort|uniq) -v -o1|cut -f2- >H.${tissue}.PA.onGene.ToRheMac8.forCmp.bed6;
			done;
		###7.1.2 Find specific PA sites based on Iso-seq reads only
			ls H.*.PA.onGene.ToRheMac8.forCmp.bed6|while read file;do
				tissue=$(echo $file|cut -f2 -d '.');
				awk -v OFS="\t" '{if($2>=30){print $1,$2-30,$3+30,$4,$5,$6}else{print $1,"0",$3+30,$4,$5,$6}}' $file|bedtools intersect -a stdin -b <(cat ../R*.PA.onGene.bed6) -v|join.pl -i1 $file -f1 4 -f2 4 -o1 |join.pl -i1 H.${tissue}.PA.onGene.bed6 -f1 4 -f2 4 >H.${tissue}SpecificPA.H.Rpos.tsv
				awk -v OFS="\t" '{if($2>=30){print $1,$2-30,$3+30,$4,$5,$6}else{print $1,"0",$3+30,$4,$5,$6}}' $file|bedtools intersect -a stdin -b ../R.${tissue}.PA.onGene.bed6 -wo|cut -f4,7-12|join.pl -i1 H.${tissue}.PA.onGene.bed6 -f1 4 -f2 1|cut -f1-6,8- |awk '{split($4,a,":");split($10,b,":");if(a[1]==b[1]){print $0}}' >H.${tissue}CommonPA.H.Rpos.tsv;
			done;
		###7.1.3 Filter by newly-defined macaque gene models
			awk -v OFS="\t" '{if($3=="+"){print $2,$5-1,$5,$12":"$5,"0",$3}else{print $2,$4,$4+1,$12":"$4+1,"0",$3}}' ~/transcriptome/visualization/Final.426909.withName.gpe |sort|uniq >Final.426909.withName.TTS.bed6 
			cat ../R*.PA.onGene.bed6 Final.426909.withName.TTS.bed6 >tmp.bed6
			ls H*SpecificPA.H.Rpos.tsv|while read file;do
				tissue=$(echo $file|cut -f2 -d '.'|sed 's/SpecificPA//');
				cut -f7-12 H.${tissue}SpecificPA.H.Rpos.tsv|perl ~/transcriptome/scripts/findHLongFromSpecific.pl -m tmp.bed6 -s >H.${tissue}Specific.longPA.bed6+ 2>H.${tissue}Specific.shortPA.bed6+
			done;
		###7.1.4 Filter using RNA-seq
			join.pl -i1 H.brainSpecific.longPA.bed6+ -f1 4 -i2 H.brainSpecific.longPA.hg19.bed6+ -f2 4 -o1|awk -v OFS="\t" '{if($6=="+"){print $1,$7,$3,$4,$5,$6}else{print $1,$2,$8,$4,$5,$6}}' >H.brainSpecific.longPA.region.rheMac8.bed6
			join.pl -i1 H.cerebellumSpecific.longPA.bed6+ -f1 4 -i2 H.cerebellumSpecific.longPA.hg19.bed6+ -f2 4 -o1|awk -v OFS="\t" '{if($6=="+"){print $1,$7,$3,$4,$5,$6}else{print $1,$2,$8,$4,$5,$6}}' >H.cerebellumSpecific.longPA.region.rheMac8.bed6
			join.pl -i1 H.heartSpecific.longPA.bed6+ -f1 4 -i2 H.heartSpecific.longPA.hg19.bed6+ -f2 4 -o1|awk -v OFS="\t" '{if($6=="+"){print $1,$7,$3,$4,$5,$6}else{print $1,$2,$8,$4,$5,$6}}' >H.heartSpecific.longPA.region.rheMac8.bed6
			ls *.region.rheMac8.bed6|sed 's/.bed6//'|while read file;do
				gpeFeature.pl -e ~/transcriptome/visualization/Final.426909.withName.gpe |bedtools subtract -a ${file}.bed6 -b stdin -s >${file}.unique.bed6
				sort -k4,4 ${file}.unique.bed6|bed6Togpe.pl|genePredToBed stdin ${file}.uniqe.bed12
			done;
			###liym@pluto:/rd1/user/liym/transcriptome/RNA-seq/rheMac8/RSeQC-FPKM/speciesSpecific/version_all/
			pythonPATH=/mnt/share/liym/python/python/Python-2.7.14
			scriptPATH=/mnt/share/liym/tools/RSeQC-2.6.3/install/home/share/local/bin
			ls *.bed6 |sed 's/.bed6//'|while read file;do bed6ToBed12.pl ${file}.bed6 >${file}.bed12; done;
			ls *bed12|while read file;do
					tissue=$(echo $file|cut -f2 -d '.'|sed 's/Specific//');
					$pythonPATH/python $scriptPATH/FPKM_count.py -i ../../../${tissue}.uniq.sorted.bam -o ${tissue} -r $file -d 1+-,1-+,2++,2-- -s 0 >log/${tissue}.log 2>log/${tissue}.err &
					$pythonPATH/python $scriptPATH/FPKM_count.py -i ../../../old.${tissue}.uniq.sorted.bam -o old.${tissue} -r $file -d 1+-,1-+,2++,2-- -s 0 >log/old.${tissue}.log 2>log/old.${tissue}.err &
					if [[ $tissue =~ "cere" ]];then
							 ls /rd1/user/lisx/annotation/RNAseq/cere-12.13/*/*/uniq.sorted.bam|while read bamFile;do
									sample=$(echo $bamFile|cut -f9 -d '/');            
									$pythonPATH/python $scriptPATH/FPKM_count.py -i $bamFile -o ${tissue}-${sample} -r $file -s 0 >log/${tissue}-${sample}.log 2>log/${tissue}-${sample}.err
							 done;
					else
							ls /rd1/user/lisx/annotation/RNAseq/${tissue}-12.13/*/*/uniq.sorted.bam|while read bamFile;do
									sample=$(echo $bamFile|cut -f9 -d '/');
									$pythonPATH/python $scriptPATH/FPKM_count.py -i $bamFile -o ${tissue}-${sample} -r $file -s 0 >log/${tissue}-${sample}.log 2>log/${tissue}-${sample}.err
							done;
					fi
			done;
			paste *brain*xls|grep -v "#"|awk '($9+$18+$27+$36+$45+$54)/6<0.2{print $4}' >brain.specificFilterbyRNA-seq.list
			paste *cerebellum*xls|grep -v "#"|awk '($9+$18+$27+$36+$45+$54)/6<0.2{print $4}' >cerebellum.specificFilterbyRNA-seq.list
			paste *heart*xls|grep -v "#"|awk '($9+$18+$27+$36+$45+$54)/6<0.2{print $4}' >heart.specificFilterbyRNA-seq.list
		###7.1.5 Final human-specific distal PA sites
			cat <(join.pl -i1 H.brainSpecific.longPA.hg19.bed6+ -f1 4 -i2 brain.specificFilterbyRNA-seq.list -o1|awk '{print $0"\tbrain"}') <(join.pl -i1 H.cerebellumSpecific.longPA.hg19.bed6+ -f1 4 -i2 cerebellum.specificFilterbyRNA-seq.list -o1|awk '{print $0"\tcerebellum"}') <(join.pl -i1 H.heartSpecific.longPA.hg19.bed6+ -f1 4 -i2 heart.specificFilterbyRNA-seq.list -o1|awk '{print $0"\theart"}') >H.specific.longPA.final.hg19.bed6+
			cat <(awk '$9=="brain"' H.specific.longPA.final.hg19.bed6+|join.pl -i1 H.brainSpecific.longPA.bed6+ -f1 4 -f2 4 -o1|awk '{print $0"\tbrain"}') <(awk '$9=="cerebellum"' H.specific.longPA.final.hg19.bed6+|join.pl -i1 H.cerebellumSpecific.longPA.bed6+ -f1 4 -f2 4 -o1|awk '{print $0"\tcerebellum"}') <(awk '$9=="heart"' H.specific.longPA.final.hg19.bed6+|join.pl -i1 H.heartSpecific.longPA.bed6+ -f1 4 -f2 4 -o1|awk '{print $0"\theart"}') >H.specific.longPA.final.rheMac8.bed6+
	##7.2 Motif comparison (Fig.5c, d)
		paste <(perl /mnt/share/liym/bin/PAmotifFinding.pl -b H.specific.longPA.final.hg19.bed6+ -l 50 -f /mnt/share/liym/data/genome/hg19/hg19.fa|cut -f1-10) <(perl /mnt/share/liym/bin/PAmotifFinding.pl -b H.specific.longPA.final.rheMac8.bed6+ -l 50 -f /mnt/share/liym/data/genome/rheMac8/rheMac8.fa|cut -f10) <(perl /mnt/share/liym/bin/PAmotifFinding.pl -b H.specific.longPA.final.hg19.bed6+ -l 50 -f /mnt/share/liym/data/ancestor/hg19_inhouse/hg19.ancestor.3species.fa|cut -f10) >H.specific.longPA.final.PAmotif.long.bed6+
		paste <(awk '{if($10!="AATAAA" && $10!="NA"){print "Other"}else{print $10}}' H.specific.longPA.final.PAmotif.long.bed6+|sort|uniq -c|awk '{print $2"\t"$1}'|sort -k1,1) <(awk '{if($11!="AATAAA" && $11!="NA"){print "Other"}else{print $11}}' H.specific.longPA.final.PAmotif.long.bed6+|sort|uniq -c|awk '{print $2"\t"$1}'|sort -k1,1|cut -f2) <(awk '{if($12!="AATAAA" && $12!="NA"){print "Other"}else{print $12}}' H.specific.longPA.final.PAmotif.long.bed6+|sort|uniq -c|awk '{print $2"\t"$1}'|sort -k1,1|cut -f2)|Rscript ~/transcriptome/scripts/barplot.PAmotif.R -names="Human,macaque,ancestor" -combine=T -o=H.specific.longPA.final.HRA.PAmotif.long.pdf
		cut -f9 H.specific.longPA.final.hg19.bed6+|sort|uniq|while read tissue;do
			awk -v OFS="\t" '{print $1,$7,$8,$4,$5,$6,$9}' H.specific.longPA.final.hg19.bed6+|awk -v value=$tissue '$7==value'|bedtools closest -s -d -a stdin -b H.${tissue}.PA.onGene.bed6 |awk -v OFS="\t" '$14<=30{split($4,a,":");split($11,b,":");if(a[1]==b[1]){print $8,$9,$10,$4,$12,$13}}'|join.pl -f1 4 -i2 <(awk -v value=$tissue -v OFS="\t" '$9==value{print $1,$7,$8,$4,$5,$6,$9}' H.specific.longPA.final.rheMac8.bed6+) -f2 4 >>H.specific.longPA.final.commonShortPA.HRlocus.tsv
		done;
		paste <(cut -f1-6 H.specific.longPA.final.commonShortPA.HRlocus.tsv|perl /mnt/share/liym/bin/PAmotifFinding.pl -l 50 -f /mnt/share/liym/data/genome/hg19/hg19.fa|cut -f1-7) <(cut -f7- H.specific.longPA.final.commonShortPA.HRlocus.tsv|perl /mnt/share/liym/bin/PAmotifFinding.pl -l 50 -f /mnt/share/liym/data/genome/rheMac8/rheMac8.fa|cut -f1-8) >H.specific.longPA.final.commonShortPA.HRmotif.tsv
		cut -f1-6 H.specific.longPA.final.commonShortPA.HRlocus.tsv|perl /mnt/share/liym/bin/PAmotifFinding.pl -l 50 -f /mnt/share/liym/data/ancestor/hg19_inhouse/hg19.ancestor.3species.fa >H.specific.longPA.final.commonShortPA.ancestor.PAmotif.bed6+
		paste <(cut -f7 H.specific.longPA.final.commonShortPA.HRmotif.tsv|awk '{if($0!="AATAAA" && $0!="NA"){print "Other"}else{print $0}}'|sort|uniq -c|awk '{print $2"\t"$1}'|sort -k1,1) <(cut -f15 H.specific.longPA.final.commonShortPA.HRmotif.tsv|awk '{if($0!="AATAAA" && $0!="NA"){print "Other"}else{print $0}}'|sort|uniq -c|awk '{print $2"\t"$1}'|sort -k1,1|cut -f2) <(cut -f7 H.specific.longPA.final.commonShortPA.ancestor.PAmotif.bed6+|awk '{if($0!="AATAAA" && $0!="NA"){print "Other"}else{print $0}}'|sort|uniq -c|awk '{print $2"\t"$1}'|sort -k1,1|cut -f2)|Rscript ~/transcriptome/scripts/barplot.PAmotif.R -names="Human,macaque,ancestor" -combine=T -o=H.specific.longPA.final.commonShortPA.HRmotif.withAncestor.pdf
	##7.3 Gene expression comparison (Fig.6a)
		###7.3.1 Mean FPKM for brain, cerebellum, heart in human and macaque
			awk -v FS="\t" -v OFS="\t" '{print $2,$16,$14,$36}' /mnt/share/liym/data/GTEx/GTEx_phe000006.v2.average_tissue_rpkm.tsv >GTEx_phe000006.v2.average_tissue_rpkm.BCH.tsv  
			####liym@pluto:/rd1/user/liym/transcriptome/RNA-seq/rheMac8/RSeQC-FPKM
			pythonPATH=/mnt/share/liym/python/python/Python-2.7.14
			scriptPATH=/mnt/share/liym/tools/RSeQC-2.6.3/install/home/share/local/bin
			$pythonPATH/python $scriptPATH/infer_experiment.py -i ../brain.uniq.sorted.bam -s 10000000 -r ../../../Final.426909.withName.bed12 >infer_experiment.brain.tsv 2>infer_experiment.brain.err
			$pythonPATH/python $scriptPATH/FPKM_count.py -i ../brain.uniq.sorted.bam -o brain -r ../../../Final.426909.withName.rmNC.merged.bed12 -d 1+-,1-+,2++,2-- -s 0 >log/brain.log 2>log/brain.err
			$pythonPATH/python $scriptPATH/FPKM_count.py -i ../cerebellum.uniq.sorted.bam -o cerebellum -r ../../../Final.426909.withName.rmNC.merged.bed12 -d 1+-,1-+,2++,2-- -s 0 >log/cerebellum.log 2>log/cerebellum.err
			$pythonPATH/python $scriptPATH/FPKM_count.py -i ../heart.uniq.sorted.bam -o heart -r ../../../Final.426909.withName.rmNC.merged.bed12 -d 1+-,1-+,2++,2-- -s 0 >log/heart.log 2>log/heart.err
			$pythonPATH/python $scriptPATH/FPKM_count.py -i ../old.brain.uniq.sorted.bam -o old.brain -r ../../../Final.426909.withName.rmNC.merged.bed12 -d 1+-,1-+,2++,2-- -s 0 >log/old.brain.log 2>log/old.brain.err
			$pythonPATH/python $scriptPATH/FPKM_count.py -i ../old.cerebellum.uniq.sorted.bam -o old.cerebellum -r ../../../Final.426909.withName.rmNC.merged.bed12 -d 1+-,1-+,2++,2-- -s 0 >log/old.cerebellum.log 2>log/old.cerebellum.err
			$pythonPATH/python $scriptPATH/FPKM_count.py -i ../old.heart.uniq.sorted.bam -o old.heart -r ../../../Final.426909.withName.rmNC.merged.bed12 -d 1+-,1-+,2++,2-- -s 0 >log/old.heart.log 2>log/heart.err
			paste *FPKM.xls|grep -v "#"|cut -f1-6,9,18,27,36,45,54|awk -v OFS="\t" '{print $4,$1":"$6":"$2"-"$3,($7+$10)/2,($8+$11)/2,($9+$12)/2}' >BCH.meanFPKM.tsv 
			ls /rd1/user/lisx/annotation/RNAseq/*-12.13/*/*/uniq.sorted.bam|while read file;do
				tissue=$(echo $file|cut -f7 -d '/'|cut -f1 -d '-');
				sample=$(echo $file|cut -f9 -d '/');
				$pythonPATH/python $scriptPATH/FPKM_count.py -i $file -o ${tissue}-${sample} -r ../../../Final.426909.withName.rmNC.merged.bed12 -s 0 >log/${tissue}-${sample}.log 2>log/${tissue}-${sample}.err
			done;
			paste *FPKM.xls|grep -v "#"|awk -v OFS="\t" '{print $4,$1":"$6":"$2"-"$3,($9+$18+$27+$36+$45+$144)/6,($54+$63+$72+$81+$90+$153)/6,($99+$108+$117+$126+$135+$162)/6}' >BCH.meanFPKM.tsv
			scp liym@pluto:/rd1/user/liym/transcriptome/RNA-seq/rheMac8/RSeQC-FPKM/BCH.meanFPKM.tsv rheMac8.newMonkey.BCH.rpkm.tsv
			perl /mnt/share/liym/bin/unique.pl -i rheMac8.newMonkey.BCH.rpkm.tsv -c 1 >rheMac8.newMonkey.BCH.rpkm.uniq.tsv
		###7.3.2 Gene expression comparison between human and macaque for genes with human-specific PA sites
			ls H*Specific.*PA.bed6+|sed 's/.bed6+//'|while read file;do
				tissue=$(echo $file|cut -f2 -d '.'|sed 's/Specific//');
				awk '$5>1' ${file}.bed6+|join.pl -f1 4 -i2 H.${tissue}.PA.onGene.bed6 -f2 4|awk -v OFS="\t" '{split($4,a,":");print a[1]":"$8,$9,$10,$11,$12,$13,$14,$7,$8}'|sort|uniq|join.pl -i2 <(cat ../R.*PA.onGene.ToHg19.bed6|sort|uniq) -f2 4 |cut -f2-7,10-|awk '$1==$7 && $6==$12'|cut -f1-6,8,9|sort|uniq >${file}.hg19.bed6+
			done;
			brain=2
			cerebellum=3
			heart=4
			rm H.commonPA.HR.FPKM.tsv H.specificLongPA.HR.FPKM.tsv H.specificShortPA.HR.FPKM.tsv
			ls H.*CommonPA.Hpos.HRmotif.bed6+|while read file;do 
				tissue=$(echo $file|cut -f2 -d '.'|sed 's/CommonPA//');
				column=$(eval echo '$'"$tissue");
				Rcolumn=$(($column+3));
				comm -23 <(awk '$5>1{print $4}' $file|cut -f1 -d ":"|sort|uniq) <(cat H.*SpecificPA.H.Rpos.tsv|cut -f4|cut -f1 -d ":"|sort|uniq)|join.pl -i1 ../GTEx_phe000006.v2.average_tissue_rpkm.BCH.tsv -o1|cut -f1,$column|join.pl -i2 ../rheMac8.newMonkey.BCH.rpkm.uniq.tsv|awk -v value=$Rcolumn -v OFS="\t" -v tissueO=$tissue '{print tissueO,$1,$2,$value}' >>H.commonPA.HR.FPKM.tsv;
				comm -23 <(awk '$5>1{print $4}' H.${tissue}Specific.longPA.hg19.bed6+|cut -f1 -d ":"|sort|uniq) <(cut -f4 H.${tissue}Specific.shortPA.hg19.bed6+|cut -f1 -d ":"|sort|uniq)|join.pl -i1 ../GTEx_phe000006.v2.average_tissue_rpkm.BCH.tsv -o1|cut -f1,$column|join.pl -i2 ../rheMac8.newMonkey.BCH.rpkm.uniq.tsv|awk -v value=$Rcolumn -v OFS="\t" -v tissueO=$tissue '{print tissueO,$1,$2,$value}' >>H.specificLongPA.HR.FPKM.tsv;
				comm -23 <(awk '$5>1{print $4}' H.${tissue}Specific.shortPA.hg19.bed6+|cut -f1 -d ":"|sort|uniq) <(cut -f4 H.${tissue}Specific.longPA.hg19.bed6+|cut -f1 -d ":"|sort|uniq)|join.pl -i1 ../GTEx_phe000006.v2.average_tissue_rpkm.BCH.tsv -o1|cut -f1,$column|join.pl -i2 ../rheMac8.newMonkey.BCH.rpkm.uniq.tsv|awk -v value=$Rcolumn -v OFS="\t" -v tissueO=$tissue '{print tissueO,$1,$2,$value}' >>H.specificShortPA.HR.FPKM.tsv;
			done;
			cat <(awk '$1=="brain"' H.specificLongPA.HR.FPKM.tsv|join.pl -f1 2 -i2 <(cut -f1 -d ":" brain.specificFilterbyRNA-seq.list) -o1) <(awk '$1=="cerebellum"' H.specificLongPA.HR.FPKM.tsv|join.pl -f1 2 -i2 <(cut -f1 -d ":" cerebellum.specificFilterbyRNA-seq.list) -o1) <(awk '$1=="heart"' H.specificLongPA.HR.FPKM.tsv|join.pl -f1 2 -i2 <(cut -f1 -d ":" heart.specificFilterbyRNA-seq.list) -o1) >H.specificLongPA.HR.FPKM.final.tsv
			cat <(awk '$1=="brain"' H.specificLongPA.HR.FPKM.PAusage.tsv|join.pl -f1 5 -i2 brain.specificFilterbyRNA-seq.list -o1) <(awk '$1=="cerebellum"' H.specificLongPA.HR.FPKM.PAusage.tsv|join.pl -f1 5 -i2 cerebellum.specificFilterbyRNA-seq.list -o1) <(awk '$1=="heart"' H.specificLongPA.HR.FPKM.PAusage.tsv|join.pl -f1 5 -i2 heart.specificFilterbyRNA-seq.list -o1) >H.specificLongPA.HR.FPKM.PAusage.final.tsv
			cut -f2 H.specificLongPA.HR.FPKM.final.tsv|sort|uniq >H.specificLongPA.HR.final.list.txt
			cut -f2 H.commonPA.HR.FPKM.tsv|sort|uniq >H.commonPA.HR.geneList.txt
	##7.4 Gene tissue distribution comparison (Fig.6b)
		awk -v FS="\t" -v OFS="\t" '{print $2,$16,$14,$36,$37,$41,$52}' /mnt/share/liym/data/GTEx/GTEx_phe000006.v2.average_tissue_rpkm.tsv >GTEx_phe000006.v2.average_tissue_rpkm.BCHKMT.tsv
		###liym@pluto:/rd1/user/liym/transcriptome/RNA-seq/rheMac8/RSeQC-FPKM
		paste BCH.meanFPKM.tsv <(paste kidney.FPKM.xls old.kidney.FPKM.xls|grep -v "#"|awk '{print ($9+$18)/2}') <(paste muscle.FPKM.xls old.muscle.FPKM.xls|grep -v "#"|awk '{print ($9+$18)/2}') <(paste testis.FPKM.xls old.testis.FPKM.xls|grep -v "#"|awk '{print ($9+$18)/2}') >BCHKMT.meanFPKM.tsv    
		scp liym@pluto:/rd1/user/liym/transcriptome/RNA-seq/rheMac8/RSeQC-FPKM/BCHKMT.meanFPKM.tsv .
		perl /mnt/share/liym/bin/unique.pl -i BCHKMT.meanFPKM.tsv -c 1 >BCHKMT.meanFPKM.uniq.tsv
		join.pl -i1 GTEx_phe000006.v2.average_tissue_rpkm.BCHKMT.tsv -i2 BCHKMT.meanFPKM.uniq.tsv >all.HR.orthologies.BCHKMT.tsv
		join.pl -i1 all.HR.orthologies.BCHKMT.tsv -i2 H.commonPA.HR.geneList.txt -o1 >HR.common.BCHKMT.FPKM.tsv
		join.pl -i1 all.HR.orthologies.BCHKMT.tsv -i2 H.specificLongPA.HR.final.list.txt -o1 >H.specificLongPA.BCHKMT.FPKM.tsv
		Rscript /mnt/share/liym/bin/rowCor.R -i1=H.specificLongPA.BCHKMT.FPKM.tsv -s1=2,7 -s2=10,15 -o=H.specificLongPA.BCHKMT.FPKM.HRcor.tsv
		paste H.specificLongPA.BCHKMT.FPKM.tsv H.specificLongPA.BCHKMT.FPKM.HRcor.tsv>tmp;mv tmp H.specificLongPA.BCHKMT.FPKM.HRcor.tsv
		Rscript /mnt/share/liym/bin/rowCor.R -i1=H.specificLongPA.BCHKMT.FPKM.tsv -s1=2,7 -s2=10,15 -m="spearman" -o=H.specificLongPA.BCHKMT.FPKM.HRSpearmanCor.tsv
		paste H.specificLongPA.BCHKMT.FPKM.tsv H.specificLongPA.BCHKMT.FPKM.HRSpearmanCor.tsv>tmp;mv tmp H.specificLongPA.BCHKMT.FPKM.HRSpearmanCor.tsv
		Rscript /mnt/share/liym/bin/rowCor.R -i1=HR.common.BCHKMT.FPKM.tsv -s1=2,7 -s2=10,15 -o=HR.common.BCHKMT.FPKM.HRcor.tsv
		Rscript /mnt/share/liym/bin/rowCor.R -i1=HR.common.BCHKMT.FPKM.tsv -s1=2,7 -s2=10,15 -m="spearman" -o=HR.common.BCHKMT.FPKM.HRSpearmanCor.tsv
		Rscript /mnt/share/liym/bin/rowCor.R -i1=all.HR.orthologies.BCHKMT.tsv -s1=2,7 -s2=10,15 -o=all.HR.orthologies.BCHKMT.HRcor.tsv
		Rscript /mnt/share/liym/bin/rowCor.R -i1=all.HR.orthologies.BCHKMT.tsv -s1=2,7 -s2=10,15 -m="spearman" -o=all.HR.orthologies.BCHKMT.HRSpearmanCor.tsv		
