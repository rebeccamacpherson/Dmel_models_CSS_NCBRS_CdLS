README.txt for MacPherson, R.A. et al., "Genetic and Genomic Analyses of Drosophila melanogaster Models of Rare Human Diseases"
Relevant GEO Accession number: GSE213763

Directory: CSS_NCBRS_CdLS_FlyModel

Files

sleep_activity_foranalysis
-- directory containing processed sleep and activity phenotype files, after elimination of dead flies, separation of day/night and activity/sleep from ShinyR-DAM output files.

sleep_activity_raw
-- directory containing raw Drosophila Acivity Monitor .txt files, trimmed for dates of experimentation, as well as a metadata file.

brains_analyses
-- directory containg raw and normalized brain dimensions as well as counts of gross abnormalities

code_behavior-brain_analyses.txt
-- file containing SAS code and R code for analysis of sleep and activity phenotypes, startle response, tapping, brain lobe morphometry, groos mushroom body abnormalities

flymodel_RNAseq_fastq_to_combinedcounts_code.txt
-- code used for merging lanes, alignment to reference genome, generation of BAM files, counting reads, and QC generation. Output file from this code is the combined_counts file in the GEO repository.

qPCR_knockdown.xlsx
-- quantitative Real Time PCR data to measure knockdown of focal genes, including delta-delta ct values for control and knockdown lines

Startle.csv
-- file containing raw startle-induced locomotor response data, with time (seconds) the fly stays moving

Tapping.csv
-- file containing a binary indication of whether the fly exhibited tapping behavior during startle reponse.

