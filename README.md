# Macaque-transcriptome
This repository contains pipeline and scripts that were used to generate resutls in paper titled "Polyadenylation-Related Isoform-switching in Human Evolution Revealed by Full-length Transcript Structure".


transcriptomeProject.sh: This file contains all the steps to analyze all the data.


/scripts: this directory contains some useful scripts for both NGS and Iso-seq data analysis. 

Reads processing:
  
	PAratio.OnGenes.sh: find PA sites based on Iso-seq reads, assign them to genes and calculate the PA usage/ratio (PA supported reads/all reads assigned to the gene) for each gene.
	dnaOrInternalPrimingContFilter.pl: remove Iso-seq reads with possible DNA or internal primer comtaminations.
	separatePaReadsFromNotFL.pl: separate Iso-seq reads with poly(A) tails from those non-full-length ones.
	fqTrimPA.pl: trim the poly(A) tail from Iso-seq reads.
	paCluster_revision.pl: output PA clusters based on Iso-seq reads in bed12 format.
	pacBioReadAssign.sh: output the pacbio reads number for each input reference transcript.
	readGroup_YM.pl: classify Iso-seq reads into unique, ambiguous or novel based on assigned gene number.
	readsAssigner.pl: assign Iso-seq reads to genes.
  
File convert or processing:
	
	geneCollaspeForRSeQC.pl: collapse all transcripts to a single transcript model for each gene based on procedures using by GTEx (https://gtexportal.org/home/documentationPage#staticTextAnalysisMethods)
	bed2gpe.pl: convert bed12 format file to genePred format file.
	sam2bed.pl: convert sam format file to bed12 format.
	filter.pl: filter a file based the input target file.
	select.pl: pick out the reference which contain the query.
	unique.pl: output unique result based on the values of specific column.
	join.pl: join two files based on sepcific columns.
	gpeFeature.pl: extract features from genePred format file.
	grep_exonfasta_from_gpe.pl: get fa format sequences for exons from genePred file.
	Base.faConvertProtein.fa.pl: convert nucleotides sequence fa format file to amino acid sequence fa format file.
  
Statistics for files:
	
	nucleotide_count_AcrossPoint.pl: count each type of nucleotides flanking the input points. 
	rowCor.R: calculate correlation value for the two input files by row.
	PAmotifFinding.pl: output 6-bp PA signal motif and its coordinate.
	PAmotifFreq.agg.pl: summarize the frequency of the input PA signals in region upstream of PA sites.
	findHLongFromSpecific.pl: output species-specific PA events.

Data visualization:

	barplot.R: generate barplots.
	barplot.PAmotif.R: generate barplot for different PA signal motif.
	hist.R: generete histograms. 
  
