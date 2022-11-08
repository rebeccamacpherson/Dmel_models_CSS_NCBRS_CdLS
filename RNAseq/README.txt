README.txt for MacPherson, R.A. et al., "Genetic and Genomic Analyses of Drosophila melanogaster Models of Chromatin Modification Disorders"
Relevant GEO Accession number: GSE213763

Directory: RNAseq

Files

RNAseq_combinedcounts-to-SASanalysis.R
-- contains R code used to normalize and filter combinedcounts file, SAS code to analyze data for differential expression, and R code to then apply FDR correction.

RNAseq_fastq-to-combinedcounts.txt
-- script used to process raw fastq files and generate a combinedcounts file (which is uploaded to GEO repository), with merging, alignment, QC, BAM file generation and read counting.

