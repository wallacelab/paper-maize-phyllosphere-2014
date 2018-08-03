#! /bin/bash

# Master script to run analysis on the 16s leaf data from the 282 in 2014 Aurora field

# Software assumed to be available via PATH (later versions are probably okay):
#   ClustalOmega 1.2.4
#   fastq-join from ea-utils (https://expressionanalysis.github.io/ea-utils/)
#   FASTX-toolkit 0.0.14
#   GNU parallel 20130922 
#   BIOM 2.1.5
#   QIIME 1.9.1
#   R 3.3.3 with packages argparse 1.0.2 , car 2.1-4, colorspace 1.2-6, corrplot 0.77, gplots 3.0.1, Hotelling 1.0-4, KEGGREST 1.12.3, lme4 1.1-12, parallel 3.4.0, and plotrix 3.6-4.
#   Python 3.4.3 with packages  argparse 1.1, biokit 0.1.2, biom-format 2.1.5, collections, gzip, itertools, json 2.0.9, math, matplotlib 1.5.0, networkx 1.11, numpy 1.10.1, os, pandas 0.16.2, re 2.2.1, scipy 0.16.1, seaborn 0.8.1 , statsmodels 0.6.1, string, and sys

# Other software
KRONA=$HOME/Software/Metagenomics/KronaTools-2.7/bin/ktImportText   # KronaTools 2.7, from https://github.com/marbl/Krona/wiki/KronaTools, down
TASSEL5="perl $HOME/Software/TASSEL/tassel-5-standalone/run_pipeline.pl -Xms10g -Xmx40g" # TASSEL verion 5.2.26
AEA=$HOME/Software/Metagenomics/AnnotationEnrichmentAnalysis/AEA    # Annotation enrichnment analysis, from Glass & Grivan 2014, DOI: 10.1038/srep04191

# Directories
rawdir=0_RawData    # Raw data
scriptdir=0_Scripts # Pipeline scripts
qualdir=$rawdir/0a_QualityCheck  # Quality control checks of raw reads
genodir=0b_MaizeGenotypes        # Maize genotypes for heritability and GWAS analysis
seqdir=1_ParseRawSequence      # Parse raw sequence into trimmed, pair-joined contigs
qiimedir=2_QiimeOtus
divdir=3_Diversity
gwasdir=4_GWAS
parsedir=5_ParseGwas
simdir=6_PowerSimulation

if [ ! -e $qualdir ]; then mkdir $qualdir; fi
if [ ! -e $genodir ]; then mkdir $genodir; fi
if [ ! -e $seqdir ]; then mkdir $seqdir; fi
if [ ! -e $qiimedir ]; then mkdir $qiimedir; fi
if [ ! -e $divdir ]; then mkdir $divdir; fi
if [ ! -e $gwasdir ]; then mkdir $gwasdir; fi
if [ ! -e $parsedir ]; then mkdir $gwasdir; fi
if [ ! -e $simdir ]; then mkdir $simdir; fi


# Global variables
maxprocs=7  # Maximum number of parallel processes
ignore_samples="LMAD_8_14ABLANK_KAK7_G5" # Remove this sample because is a known blank that failed. (Had reads, but didn't look like any of the plant samples.) 
plate_key=$scriptdir/0_plate_key.txt
chromlengths=$scriptdir/0_Zmays_agpv3_chrom_lengths.txt

# # # TODO : delete? transcript_data=0_KarlData/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_rounded.txt
# # # TODO: delete? transcript_key=0_KarlData/laptop_version_Master_name_conversion_from_orig_1960_first_sheet_for_Wallace_01042016.txt


# # # Do quality checking on the HiSeq reads
# # NOTE: 1 sample will throw a mismatched name error (LMAD_26_14ABLANK_KAK12_G1) because the server cut out the first part of the sample name for some reason
# md5_sums=$scriptdir/0_md5_sums.txt
# ./0a_QualityCheckHiseqReads.sh $qualdir $rawdir $scriptdir $plate_key $md5_sums

# # Gather the maize genotype data for later heritability and GWAS analysis
# # NOTE: This is mostly included for reference; the filtered genotype data should be included as part of the pipeline download
# source_genos=$HOME/Projects/0_RawData/MaizeGBS/Panzea_v2.7/AllZeaGBSv2.7_publicSamples_imputedV3b_agpv3.hmp.gz
# source_keyfile=$HOME/Projects/0_RawData/MaizeGBS/Panzea_v2.7/AllZeaGBSv2.7_publicSamples_metadata20140411.txt
# additional_samples=$scriptdir/0_additional_samples_manual.txt
# synonyms="b2-good:b2 ci28a:ci28agoodman-buckler ci.7:ci7goodman-buckler ci91b:ci91bgoodman-buckler de_2:de2 il101:il101t hy:ill.hy mog:mo.g mr19_(santo_domingo):mr19santodomingo mr20_(shoe_peg):mr20shoepeg w22r-r-std_cs-2909-1:w22r-rstd"
# min_freq=0.05   # Filter out sites with minor allele frequency below this
# min_count=100   # filter out sites that don't have genotypes for at least this many samples
# ./0b_GatherMaizeGenotypes.sh $scriptdir $genodir $plate_key $source_genos $source_keyfile $additional_samples "$synonyms" "$TASSEL5" $min_freq $min_count

# # # Parse raw reads - quality filter, join paired ends, and add sample names
# ./1_ParseRawReads.sh $rawdir $scriptdir $seqdir 

# # QIIME OTU picking
# min_sample_abundance_over_blanks=1.1  # Ratio of a sample's reads to the highest blank for it to pass filtering
# read_dist=$seqdir/1b_read_counts.txt
# min_otu_size=10    # Minimum number of times an OTU has to be seen to be picked by QIIME (default is 2)
# ./2_PickQiimeOtus.sh $scriptdir $qiimedir $seqdir $plate_key $read_dist $min_otu_size $min_sample_abundance_over_blanks $ignore_samples $maxprocs $KRONA $rawdir

# # # Comunity diversity statistics; also determines the core community
# prefix=$qiimedir/2f_otu_table.sample_filtered.no_mitochondria_chloroplast
# rarefaction_level=5000
# core_fraction=0.8   # Fraction of samples an OTU has to be in for it to count as core
# ./3_DiversityStatistics.sh $scriptdir $divdir $prefix $rarefaction_level $core_fraction $KRONA $maxprocs


# # # TODO: Check beta-diversity with top 30 or top 50 OTUs and see if is any more heritable?
# # Calculate heritability and run GWAS on highly heritable OTUs, disversity statistics, and implied metagenome content
# qiime_prefix=$qiimedir/2f_otu_table.sample_filtered.no_mitochondria_chloroplast
# rarefaction_level=10000   # Minimum reads for a sample to be included
# min_counts_across_samples=0.5  # How many samples an OTU (or other unit) has to have at least 1 read in.
# num_pcs=5   # Number of principal coordinates to include in analysis
# flowering_time=$scriptdir/0_flowering_times.txt
# cofactors="tissue_date collector"
# genotypes=$genodir/0h_samples_sorted_filtered.hmp.txt.gz
# aliasfile=$scriptdir/0_inbred_aliases_manual.txt
# top_n=500   # Focus on top N OTUs
# nperms=100  # Number of permutations for heritability analysis (first pass; p-vlaue of 0.05 used as filter)
# bigperms=10000 # Number of permutations to the ones that pass initial fitlering from above line
# perm_cutoff=0.001 # Significance threshold for the bigperms empirical p-value (initial one always at 0.02)
# gwas_perms=100  # number of permutations to determine significance of GWAS
# gwas_procs=4    # Number of processor cores to use for GWAS (since takes more memory than most other processes, best to have this smaller)
# cogdir=$scriptdir/COG_data  # Directory holding COG annotation files for metagenome annotation
# empirical_pval_cutoff=0.1 # Empirical p-value cutoff for permutation analysis (usually err on the side of generosity, so ~0.1)
# all_flowering=$scriptdir/0_goodman_average_flowering.panzea.txt # Average flowering time used to interpolate flowering values to get more power in GWAS
# cluster_cutoff=0.3  # Linkage disequilibrium cutoff to consider SNP hits in different linkage clusters (to remove redundant SNPs)
# ./4_RunHeritabilityAndGwas.sh $scriptdir $gwasdir $qiime_prefix $rarefaction_level $min_counts_across_samples $num_pcs $flowering_time "$cofactors" $genotypes $aliasfile $top_n $nperms $bigperms $perm_cutoff \
#   $gwas_perms $chromlengths "$TASSEL5" $cogdir $maxprocs $gwas_procs $empirical_pval_cutoff $all_flowering $cluster_cutoff


# # Parse the GWAS results
# qiime_prefix=$qiimedir/2f_otu_table.sample_filtered.no_mitochondria_chloroplast
# herit_pval_cutoff=0.001
# max_h2=0.98 # Some traits have really, really high heritabilities that are probably artifacts, so filter out
# h2_perms=100000  # number of random permutations for metagenome Annotation Enrichment Analysis
# blups=$gwasdir/4u_combined_good_traits.update_flowering.tassel.txt  # file of BLUPs to use for interpreting GWAS results
# min_cluster_score=1.5   # Number of traits (after correlation correction) that hit in the same window to look at clustering
# cluster_p=0.01  # Empirical p-value cutoff to use when analyzing clusters of hits
# ./5_ParseGwasResults.sh $scriptdir $parsedir $gwasdir $qiime_prefix $AEA $herit_pval_cutoff $max_h2 $h2_perms $chromlengths $blups $min_cluster_score $cluster_p


# Do a simulation study as requested by reviewers. (Added after initial manuscript submission)
#   NOTE: This step generally used more recent software and package versions than the above versions because a system update had taken place between submission and revisions
genotypes=$gwasdir/4u_gwas/all_genos.hmp.txt.gz
gwas_outdir=$gwasdir/4u_gwas/4u_orig_gwas/  # Directory where original GWAS results are kept
genos_bychrom_prefix=$gwasdir/4u_gwas/chr
phenotypes=$gwasdir/4u_combined_good_traits.update_flowering.tassel.txt
covariate=$gwasdir/4u_flowering_covariate.interpolated.tassel.txt
permdir=$gwasdir/4u_gwas/4u_permuted_genos
traits="flowering_time log_COG4679_3080 weighted_unifrac_PC1 log_qiime.Bacteria.Proteobacteria.Alphaproteobacteria.Rhizobiales.Methylobacteriaceae.Methylobacterium.adhaesivum.591699"   # Traits to run simulations on
# traits="log_COG4679_3080 weighted_unifrac_PC1 log_qiime.Bacteria.Proteobacteria.Alphaproteobacteria.Rhizobiales.Methylobacteriaceae.Methylobacterium.adhaesivum.591699"   # Traits to run simulations on
num_qtn="5 20 100"  # Comma-separated list of number of QTN to simulate in each dataset
herits="0.2 0.5 0.8"    # Comma-separated list of heritabilities to simulate. NOTE: This is equivalent to the fraction variance explained _after_ fitting a kinship matrix.
nreps=10   # number of replicate simulations at each QTN/heritability level
empirical_pval_cutoff=0.1
./6_DoSimulationAnalysis.sh $scriptdir $simdir $genotypes $phenotypes $covariate "$traits" "$num_qtn" "$herits" $nreps $gwas_outdir $genos_bychrom_prefix $chromlengths $permdir $empirical_pval_cutoff $maxprocs "$TASSEL5"
