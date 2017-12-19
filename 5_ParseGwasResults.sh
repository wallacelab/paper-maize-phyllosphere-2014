#! /bin/bash 

# Collapse OTU counts to BLUPs and then run statistics off of those

scriptdir=$1
parsedir=$2
gwasdir=$3
qiime_prefix=$4
AEA=$5  # annotation enrichment analysis software
herit_pval_cutoff=$6 # p-value cutoff for heritability
max_h2=$7
h2_perms=$8
chromlengths=$9
blups=${10}
min_cluster_score=${11}
cluster_p=${12}

biom=$qiime_prefix.biom
keyfile=$qiime_prefix.key.tsv
taxonomy=$qiime_prefix.taxonomy.txt
tree=$qiime_prefix.tre



# # Parse metagenome heritability to easier files
# python3 $scriptdir/5_AssembleMetagenomeHeritData.py -i $gwasdir/4s_narrow_herit_bigperms_metagenome*/*/4r_*_stats.txt -o $parsedir/5_metagenome_data #--debug
# pigz $parsedir/5_metagenome_data.perms.txt

##########
# Fisher Exact and Annotation enrichment analysis for heritability
##########

## Annotation Enrichment files:
# Index file:      3 columns, with gene id #, gene name, and total count of terms
# Annotation file: 2 columns, with GO term name and comma-separated list of gene IDs annotated to it
# Partition file:  2 columns, with GO term name and comma-separated list of other GO terms it is a child of. (First one is itself and gives its own ID.)

# for set in cog ko; do
#   prefix=$parsedir/5b_annotations.$set
  # Set up graph of annotations
#   python3 $scriptdir/5a_AssembleMetagenomeHierarchy.py -i $gwasdir/4j_predicted_metagenome.$set.biom -o $parsedir/5a_metagenome_key.$set
# # 
# #   # Format for Annotation Enrichment Analysis
#   python3 $scriptdir/5b_FormatHierarchyForAnnotationEnrichment.py -i $parsedir/5a_metagenome_key.$set.hierarchy.gexf -o $prefix --set $set
#   python3 $scriptdir/5b_FormatMetagenomeHeritForAnnotationEnrichment.py -i $parsedir/5_metagenome_data.herits.txt -s $set -o $prefix.signature.txt \
#     --pval-cutoff $herit_pval_cutoff --max-herit $max_h2  --name ${set}_enrichment
#   
# #   # Run annotation enrichment and format back into more readable
#   $AEA -i $prefix.index.txt -a $prefix.annotations.txt -p $prefix.partitions.txt -s $prefix.signature.txt -o $parsedir/5c_annotations.$set.enrichment.txt -n $h2_perms
#   python3 $scriptdir/5c_ReformatAeaOutput.py -i $parsedir/5c_annotations.$set.enrichment.txt -k $parsedir/5b_annotations.$set.term_key.txt -o $parsedir/5c_annotations.$set.enrichment.reformatted.txt
# 
#   # Now run Fisher's Exact test to compare
#   python3 $scriptdir/5d_FisherExactForMetagenome.py -i $parsedir/5_metagenome_data.herits.txt -s $set --pval-cutoff $herit_pval_cutoff --max-herit $max_h2 \
#     -o $parsedir/5d_fisher_exact.$set --graphfile $parsedir/5a_metagenome_key.$set.hierarchy.gexf

#   # Make subgraph of included terms
#   python3 $scriptdir/5e_MakeSubgraphOfTerms.py -i $parsedir/5_metagenome_data.herits.txt  -g $parsedir/5a_metagenome_key.$set.hierarchy.gexf -o $parsedir/5e_fisher_exact.$set.hits_graph \
#     --pval-cutoff $herit_pval_cutoff --max-herit $max_h2
#   break
# done

# Rscript -e "a=read.delim('$parsedir/5a_metagenome_key.cog.txt'); b=read.delim('$parsedir/5a_metagenome_key.ko.txt')" \
#   -e "write.table(rbind(a,b), file='$parsedir/5a_metagenome_key.txt', sep='\t', quote=F, row.names=F, col.names=T)"


###########
# Parse GWAS results
###########

all_hits=$gwasdir/4x_gwas_hits.clustered.all.txt
best_hits=$gwasdir/4x_gwas_hits.clustered.best.txt

# Get correlations among the BLUPs so can weight common hits
Rscript -e "x=read.delim('$blups', skip=2, row.names=1); cors=cor(x, use='pairwise');" \
  -e "write.table(cors, file='$parsedir/5f_blup_correlations.txt', sep='\t', quote=F, row.names=T, col.names=T)"

# Look for common hits in specific windows, using different p-value cutoffs
for winsize in 1000000 100000 10000; do
  # Compile hits and graph clusters
#   python3 $scriptdir/5g_PlotCommonSnpClusters.py -i $best_hits -o $parsedir/5g_common_hits.winsize$winsize --winsize $winsize --chromlengths $chromlengths --cors $parsedir/5f_blup_correlations.txt #--debug

  # Look at correlations/clusters among the multi-hit regions to see what types of traits are clustering there
  for p in 0.05 0.01; do
    python3 $scriptdir/5h_AnalyzeClusteredTraits.py -i $parsedir/5g_common_hits.winsize$winsize.hitcounts.txt -o $parsedir/5h_common_hits.winsize$winsize --blups $blups \
      --pval-cutoff $p --min-score $min_cluster_score --key $parsedir/5a_metagenome_key.txt
  done
#   break
done



# # # # Compile empirical p-value results into a single file and plot. Note that in each window, each trait only provides 1 hit since SNPs are likely in LD
# # # gwas_results=$gwasdir/4u_gwas
# # # lasso_results=$gwasdir/4w_gwas_lasso/4x_final_models
# # # 
# # # # for winsize in 1000000 100000 10000; do
# # # #   # Compile hits and graph clusters
# # # #   python3 $scriptdir/5g_CompileEmpiricalPvals.py -i $gwas_results/4v_gwas_results/*.p_adjusted.txt -o $parsedir/5g_compiled_gwas.empirical_pvals.winsize$winsize \
# # # #   --winsize $winsize --chromlengths $chromlengths #--debug
# # # # 
# # # #   # Select results from the above and turn into a table for identifying hit regions
# # # #   python3 $scriptdir/5h_PullOutHitClusters.py -i $parsedir/5g_compiled_gwas.empirical_pvals.winsize$winsize.hitcounts.txt -f $parsedir/5g_compiled_gwas.empirical_pvals.winsize$winsize.txt \
# # # #     --cutoff 0.05 --min-hits 5 -o $parsedir/5h_hit_clusters.winsize$winsize.txt -k $parsedir/5a_metagenome_key.txt
# # # # 
# # # #   # Look for correlations among the phenotypes to test if they're all just different names for the same things
# # # #   Rscript $scriptdir/5i_LookForCorrelationsInClusters.r -i $parsedir/5h_hit_clusters.winsize$winsize.txt -b "$gwas_results/4u_blups/4u_blups.%%%.txt" -o $parsedir/5i_cluster_correlation.$winsize.png
# # # 
# # # #   break
# # # # done
# # # 
# # # # Compile LASSO and empirical p-value results results
# # # 
# # # # Rscript $scriptdir/5j_CompileLassoResults.r -i $lasso_results/*.lasso_model.txt -o $parsedir/5j_lasso_models.txt --name-pattern ".+/4x_(.+).lasso_model.txt"
# # # Rscript $scriptdir/5j_CompileEmpiricalPvalResults.r -i $gwas_results/4v_gwas_results/*.p_adjusted.txt -o $parsedir/5j_empirical_pval_models.txt 
# # # 
# # # # TODO: Go back and use LASSO with functionally independent tests instead of FDR? Multiple testing is just killing me, but tests are not independent
# # # 
# # # # TODO: Runing GWAS on residuals doesn't work. Misses vgt1, whereas MLM on it hits. Unfortunately, is also much, much slower. Figure out a compromise?
# # # 
# # # # Check for correlations among traits
# # # picrust_data=/usr/local/lib/python2.7/dist-packages/PICRUSt-1.1.0-py2.7.egg/picrust/data/
# # # source_biom=$gwasdir/4i_closed_reference_otus.fix_copy_number.biom
# # # winsize=10000
# # # # biom convert -i $source_biom -o $parsedir/5k_closed_reference_otus.fix_copy_number.biom.txt --to-tsv
# # # # cut -f1 $parsedir/5k_closed_reference_otus.fix_copy_number.biom.txt | grep -v "^#" > $parsedir/5k_target_otus.txt
# # # # Rscript -e "x=read.delim('$parsedir/5h_hit_clusters.winsize$winsize.txt', stringsAsFactors=F); write(unique(x\$Trait), file='$parsedir/5k_target_traits.txt')"
# # # # python3 $scriptdir/5k_SubsetPicrustData.py -i $picrust_data/ko_13_5_precalculated.tab.gz $picrust_data/cog_13_5_precalculated.tab.gz --otus  $parsedir/5k_target_otus.txt --traits $parsedir/5k_target_traits.txt \
# # # #   -o $parsedir/5k_trimmed_picrust_data.txt #--debug
# # # 
# # # # # Look at correlations among traits just in OTUs
# # # # Rscript $scriptdir/5l_PlotCorrelationsAmongTraits_rawdata.r -i $parsedir/5k_trimmed_picrust_data.txt -o $parsedir/5l_trimmed_picrust_data.cors --hits $parsedir/5h_hit_clusters.winsize$winsize.txt \
# # # #   -s $parsedir/5k_closed_reference_otus.fix_copy_number.biom.txt
# # # 
# # # 
# # # # Step 1: Filter raw data files for just the OTUs and KO/COG terms I'm interested in to make them much smaller
# # # # Step 2: Look for correlations among terms at level of reference OTUs, and then after multiplying the matrices but not summing across samples. 
# # # #   (Basically, are these just highly correlated within organisms, or did I get a fluke where a lot are correlated after running analysis but not in the raw data)
# # # 
# # # # NOTE: So far, things just seem to be highly correlated. I've got high correlatoins (+ or -) on the marginal sums of the matrices. Still need to check the raw matrices, although the answer is almost certainly yes








# # # # Options for meta-analysis
# # # #   PLINK - Does it on its own, but would need to convert my output to PLINK's format: http://zzz.bwh.harvard.edu/plink/metaanal.shtml
# # # 
# # # # Take top 1% or 0.1% in all studies and rank by how often they're hit; bin into megabase windows? (and compare to null distribution)
# # # #   Suffers from going for specific peaks over regions. Would really want to take model-selected SNPs with enter and exit values
# # # 
# # # # Papers to check out
# # # #  http://www.nature.com/nrg/journal/v14/n6/full/nrg3472.html#ref83 - Review of metaanalysis methods. Mostly for same-pehnotype analysis, but some notes on diff-phenotype ones
# # # #  https://www.ncbi.nlm.nih.gov/pubmed/18976227 - Multivariate meta-analysis. See if can adapt
# # # #  http://onlinelibrary.wiley.com/doi/10.1002/gepi.20556/full - Two-step LASSO; maybe a form of model selection to use instead of empirical p-value? I'm liking this option. Take existing results, get FDR, do LASSO on FDR>=0.05
# # # #  http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002254#s4 - Cross-phenotype metanalysis of P-values based on them supposed to have a uniform distribution if there are no additional associations
# # # #    Potential problem is that it was done on only ~150 SNPs; would need an initial screen first b/c scaling to 400k would probably crash the computer. Highly correlated SNPs dealt with by clustering post-hoc to find
# # # #    what the correlations underlie
# # # 
# # # # cvfit <- glmnet::cv.glmnet(x, y)
# # # # coef(cvfit, s = "lambda.1se")
# # # # Maybe try multiresponse gaussian family to look at response across all; probably have too many, though.





