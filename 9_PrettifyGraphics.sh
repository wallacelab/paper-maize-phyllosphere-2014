#! /bin/bash

# This script is not technically part of the analysis pipeline. Instead, it takes the results from the pipeline and formats them into publication-ready graphics, extracts data for tables, etc.

# Software
KRONA=$HOME/Software/Metagenomics/KronaTools-2.7/bin/ktImportText

# Directories
rawdir=0_RawData    # Raw data
scriptdir=0_Scripts # Pipeline scripts
qualdir=$rawdir/0a_QualityCheck  # Quality control checks of raw reads
genodir=0b_MaizeGenotypes        # Maize genotypes for heritability and GWAS analysis
parsedir=1_ParseRawSequence      # Parse raw sequence into trimmed, pair-joined contigs
qiimedir=2_QiimeOtus
divdir=3_Diversity
gwasdir=4_GWAS
parsedir=5_ParseGwas/
plotdir=9_PrettyGraphics

# Specific subdirectories for intermediate files
kronadir=$plotdir/9a_KronaPlot
coredir=$plotdir/9c_CoreOtus
pairdir=$plotdir/9d_PairedDayNightSamples
alphadir=$plotdir/9f_AlphaDiversity
betadir=$plotdir/9h_BetaDiversity
heritdir=$plotdir/9j_Heritability

if [ ! -e $plotdir ]; then mkdir $plotdir; fi
if [ ! -e $kronadir ]; then mkdir $kronadir; fi
if [ ! -e $coredir ]; then mkdir $coredir; fi
if [ ! -e $pairdir ]; then mkdir $pairdir; fi
if [ ! -e $alphadir ]; then mkdir $alphadir; fi
if [ ! -e $betadir ]; then mkdir $betadir; fi
if [ ! -e $heritdir ]; then mkdir $heritdir; fi

# Frequently used variables
qiime_prefix=$qiimedir/2f_otu_table.sample_filtered.no_mitochondria_chloroplast # Raw, non-normalized qiime file
normalized_otus=$gwasdir/4f_orig_otus.filtered.txt # Normalized qiime file


#########################
# Krona plot of maize diversity
#########################


# # Make Krona plots of the unnormalized biom file
# python3 $scriptdir/2g_FormatTaxonomyForKrona.py -i $qiime_prefix.taxonomy.txt -b $qiime_prefix.biom  -o $kronadir/9a_krona_input.raw.txt --include-otus
# $KRONA -o $kronadir/9b_krona_plot.html $kronadir/9a_krona_input.raw.txt
# 
# # Do the same with the normalized one and compare
# biom convert -i $normalized_otus -o $kronadir/9a_normalized_otus.biom --to-hdf5
# python3 $scriptdir/2g_FormatTaxonomyForKrona.py -i $qiime_prefix.taxonomy.txt -b $kronadir/9a_normalized_otus.biom  -o $kronadir/9a_krona_input.normalized.txt --include-otus
# $KRONA -o $kronadir/9b_krona_plot.normalized.html $kronadir/9a_krona_input.normalized.txt


# At this point the "snapshot" feature was used to take an SVG snapshot of the normalized data, which was modified in Inkscape for the final image. (Sorry, Krona doesn't have a good API to tweak things.)


#########################
# Core OTUs
#########################

# # Defining the core microbiome as being ones that are present with at least 10 out of the 10,000 reads the data was normalized to, and in 80% of samples
# Rscript $scriptdir/9c_GetCoreOtuStats.r -i $normalized_otus -o $coredir/9c_core_otus.txt -a $coredir/9c_all_otus.txt --min-count-in-samples 10 --min-presence-across-samples 0.8
# 
# # Distribution of core OTUs within samples, and distribution of major bacterial divisions across samples
# python3 $scriptdir/9c_PlotOtuDistributionsWithinAndAcrossSamples.py -i $normalized_otus -t $qiime_prefix.taxonomy.txt -o $coredir/9c_otu_distributions_within_across_samples -c $scriptdir/9c_taxonomy_colorkey.txt \
#   --core-otus $coredir/9c_core_otus.txt --sample-filter "LMAD.8."


#########################
# Compare Day/Night samples of same leaves
# Only difference is that were sampled 12 hours apart on opposite sides of the leaves
#########################

# # # NOTE: All of following samples were on Aug 8
# python3 $scriptdir/9d_PlotDayNightDifferences_publication.py -i $normalized_otus -c $coredir/9c_core_otus.txt -o $pairdir/9d_core_otus_day_versus_night.publication \
#   -k $qiime_prefix.key.tsv --plot-order --lines Oh43 Ia5125 C123 A641 A679 --name-key $scriptdir/9d_otu_name_key.txt #--debug 

# bray_curtis=$divdir/3a_diversity_bray_curtis/bdiv_even5000/bray_curtis_dm.txt
# weighted=$divdir/3a_diversity_unifrac/bdiv_even5000/weighted_unifrac_dm.txt
# unweighted=$divdir/3a_diversity_unifrac/bdiv_even5000/unweighted_unifrac_dm.txt
# 
# for distances in $bray_curtis $weighted $unweighted; do
#   set=`basename $distances`
#   set=${set/_dm.txt/}
#   python3 $scriptdir/9e_PLotDistancesBetweenDayAndNightSamples.py -i $distances -o $pairdir/9e_distance_comparison.$set
# #   break
# done


#########################
# Alpha diversity
#########################

# # # Alpha diversity analysis on normalized data, and a typical rarefaction of the original data (for comparison)
# metrics=$(alpha_diversity.py --show_metrics | grep "^Known" | sed -e "s/Known metrics are: //" -e "s/ //g" -e "s/michaelis_menten_fit,//" -e "s/kempton_taylor_q,//"  )  # Michaelis-Menten hangs analysis, so skip. Kempton_taylor has issues with all NAs
# single_rarefaction.py -i $qiime_prefix.biom -o $alphadir/9f_alpha_diversity_rarefied.biom -d 10000
# alpha_diversity.py -i $alphadir/9f_alpha_diversity_rarefied.biom -o $alphadir/9f_alpha_diversity_rarefied_data.txt -m $metrics -t $qiime_prefix.tre
# alpha_diversity.py -i $normalized_otus -o $alphadir/9f_alpha_diversity_normalized_data.txt -m $metrics -t $qiime_prefix.tre

# # Prune to metrics with <0.9 correlation and then find Type 2 SS significance for cofactors
# cofactors="date time rna_plate collector subpop"
# cofactors_sub="date time rna_plate subpop"  # See what effect of Collector is, since explains so much variance by itself
# for set in normalized rarefied; do
#   Rscript $scriptdir/9g_PruneAlphaDiversityMetrics.r -i $alphadir/9f_alpha_diversity_${set}_data.txt -o $alphadir/9g_alpha_diversity_${set}_data --correlation-anchor berger_parker_d --prune 0.9
#   Rscript $scriptdir/9h_AnalyzeAlphaDiversityCofactors.r -i $alphadir/9g_alpha_diversity_${set}_data.non_redundant_metrics.txt -o $alphadir/9h_alpha_diversity_${set} -f $cofactors -m $qiime_prefix.key.tsv --hidden-factors
#   Rscript $scriptdir/9h_AnalyzeAlphaDiversityCofactors.r -i $alphadir/9g_alpha_diversity_${set}_data.non_redundant_metrics.txt -o $alphadir/9h_alpha_diversity_${set}_no_collector -f $cofactors_sub -m $qiime_prefix.key.tsv --hidden-factors
# done


# # # Based on above anlaysis, chose 13 hidden factors and include collector as a cofactor because it does have a significant effect
# # # This step is basically outputting the variance explaiend output in pretty form. Also used the variance explained plot from the previous step
# cofactors="date time rna_plate collector subpop"
# Rscript $scriptdir/9h_AnalyzeAlphaDiversityCofactors_publication.r -i $alphadir/9g_alpha_diversity_normalized_data.non_redundant_metrics.txt -o $alphadir/9i_alpha_diversity_normalized_data.pretty \
#   -f $cofactors -m $qiime_prefix.key.tsv -n 13

# rarefaction_level=9500  # Smallest value in my normalized dataset is 9664, so get close to that
# core_diversity_analyses.py -i $qiime_prefix.biom -o $alphadir/9j_qiime_diversity_orig_data -m $qiime_prefix.key.tsv --sampling_depth $rarefaction_level --tree_fp $qiime_prefix.tre --suppress_beta_diversity --suppress_taxa_summary
# core_diversity_analyses.py -i $normalized_otus -o $alphadir/9j_qiime_diversity_normalized_data -m $qiime_prefix.key.tsv --sampling_depth $rarefaction_level --tree_fp $qiime_prefix.tre --suppress_beta_diversity --suppress_taxa_summary

# Make alpha diversity plots showing increased diversity with date, subpop, and collector. 
# # targetdir=$alphadir/9j_qiime_diversity_normalized_data/arare_max$rarefaction_level
# make_rarefaction_plots.py -i $targetdir/alpha_div_collated/ -m $qiime_prefix.key.tsv -o $targetdir/alpha_rarefaction_plots --colorby date,subpop,collector --generate_average_tables
# python3 $scriptdir/9j_PlotAlphaDiversity.py -i $targetdir/alpha_rarefaction_plots/average_tables/*.txt -o $alphadir/9j_alpha_diversity.png

# TODO: test for statistical significance of alpha-diversity differences among groups?

#########################
# Beta diversity
#########################

rarefaction_level=9500  # Smallest value in my normalized dataset is 9664, so get close to that
metrics="bray_curtis weighted_unifrac unweighted_unifrac"

# # Modified pipeline from core_diversity_analyses.py from QIIME, taken out to have a bit more control over parameters and file names
# for set in normalized orig; do
#   workdir=$betadir/9g_qiime_diversity_${set}
#   if [ ! -e $workdir ] ; then mkdir $workdir; fi
#   if [ $set == "normalized" ] ; then 
#     infile=$normalized_otus 
#   else
#     infile=$qiime_prefix.biom
#   fi
#   
#   # Filter and rarefy
#   filter_samples_from_otu_table.py -i $infile -o $workdir/9f_table_mc$rarefaction_level.biom -n $rarefaction_level
#   single_rarefaction.py -i $workdir/9f_table_mc$rarefaction_level.biom -o $workdir/9f_table_even$rarefaction_level.biom -d  $rarefaction_level
#   
#   # Distance metrics and PCoA calculation
#   mymetrics=`echo $metrics | tr ' ' ','`
#   beta_diversity.py -i $workdir/9f_table_even$rarefaction_level.biom -o $workdir --metrics $mymetrics -t $qiime_prefix.tre
#   
#   # Process each metric in turn
#   for metric in $metrics; do
#     mv $workdir/${metric}_9f_table_even$rarefaction_level.txt $workdir/${metric}_dm.txt
#     principal_coordinates.py -i $workdir/${metric}_dm.txt -o $workdir/${metric}_pc.txt 
#     make_emperor.py -i $workdir/${metric}_pc.txt -o $workdir/${metric}_emperor_pcoa_plot/ -m $qiime_prefix.key.tsv 
#   done
# done

# # Some checks on the different metrics
# for metric in $metrics; do
#   # Get correlation of first 10 PCs in whole-dataset versus normalized; mostly for sanity checking
#   Rscript $scriptdir/9g_CorrelateBetaPcs.r --norm $betadir/9g_qiime_diversity_normalized/${metric}_pc.txt --orig $betadir/9g_qiime_diversity_orig/${metric}_pc.txt \
#     -o $betadir/9g_correlate_orig_norm.$metric.png -n 10 --name $metric
# 
#   #Calculate statistical significance of categories for PCA clustering
#   Rscript $scriptdir/9h_GetSignificanceOfFactorsForPcs.r -i $betadir/9g_qiime_diversity_normalized/${metric}_pc.txt -o $betadir/9h_manova.$metric -n 10 \
#     --name $metric --include-batches --keyfile $qiime_prefix.key.tsv --skip Description row col plate_row plate_col
# 
#   break
# done


# # Plot beta diversity publication graphics show separation by dna_plate, and tissue-date
# python3 $scriptdir/9h_PlotBetaDiversity.py -i $betadir/9g_qiime_diversity_normalized/*_pc.txt -o $betadir/9h_beta_diversity -k $qiime_prefix.key.tsv --categories dna_plate tissue_date



##################
# Heritability
##################

# # Gather heritability stats for the manuscript; mostly for convenience and to record which ones I'm drawing from.
# # Beta diversity
# cp $gwasdir/4q_narrow_herit_diversity.txt $heritdir/9j_diversity_herit.small.txt
# cp $gwasdir/4s_narrow_herit_bigperms_diversity.txt $heritdir/9j_diversity_herit.bigperms.txt

# # OTUs
# cp $gwasdir/4q_narrow_herit_otus.txt $heritdir/9j_otus_herit.small.txt
# cp $gwasdir/4s_narrow_herit_bigperms_otus.txt $heritdir/9j_otus_herit.bigperms.txt
# python3 $scriptdir/9j_CombineHeritTables.py --small $heritdir/9j_otus_herit.small.txt --big $heritdir/9j_otus_herit.bigperms.txt -o $heritdir/9j_otus_herit.combined.txt --type otu

# # Metagenome
metadir=$heritdir/9j_metagenome_tables
# if [ ! -e $metadir ]; then mkdir $metadir; fi
# cp $gwasdir/4q_narrow_herit_metagenome?.txt $metadir
# cp $gwasdir/4s_narrow_herit_bigperms_metagenome?.txt $metadir
# for set in ko cog; do
#   tail -q -n+2  $gwasdir/4j_predicted_metagenome.$set*.biom.txt | cut -f1 | sort | uniq > $heritdir/9j_${set}.traits.txt
#   echo "flowering_time" >> $heritdir/9j_${set}.traits.txt
#   python3 $scriptdir/9j_CombineHeritTables.py --small $metadir/4q_narrow_herit_metagenome?.txt --big $metadir/4s_narrow_herit_bigperms_metagenome?.txt \
#     -o $heritdir/9j_metagenome_herit.${set}.txt --type metagenome --metagenome-key $parsedir/5a_metagenome_key.${set}.txt --traits $heritdir/9j_${set}.traits.txt
# #   break
# done
# combine metagenome keys
# Rscript -e "ko=read.delim('$parsedir/5a_metagenome_key.ko.txt'); cog=read.delim('$parsedir/5a_metagenome_key.cog.txt'); out=rbind(ko, cog)" \
#   -e "write.table(out, file='$heritdir/9j_metagenome_key.combined.txt', sep='\t', quote=F, row.names=F, col.names=T)"
# Combine metagenome heritability tables
# python3 $scriptdir/9j_CombineHeritTables.py --small $metadir/4q_narrow_herit_metagenome?.txt --big $metadir/4s_narrow_herit_bigperms_metagenome?.txt \
#   -o $heritdir/9j_metagenome_herit.combined.txt --type metagenome --metagenome-key $heritdir/9j_metagenome_key.combined.txt --traits $heritdir/9j_*.traits.txt

# Plot heritabilities 
# cutoff=0.001
# # OTUs
# python3 $scriptdir/9k_PlotHeritability_publication.py -i $gwasdir/4s_narrow_herit_bigperms_otus/*/4r_*_stats.txt -g $heritdir/9k_otu_heritability_pretty \
#   --pval-cutoff $cutoff -n 20 --figsize 5 6 #--debug

# Separate heritability plotting into 2 steps so don't have to load everything each time
# python3 $scriptdir/5_AssembleMetagenomeHeritData.py -i $gwasdir/4s_narrow_herit_bigperms_metagenome*/*/4r_*_stats.txt -o $heritdir/9k_metagenome_data #--debug


###########
# Metagenome heritability graphs, etc.
###########

# Copy over Fisher Exact and Annotation Enrichment files
# cutoff=0.001
# max_h2=0.98 # Things above this are thought to be artifacts due to their distribution os permuted heritabilities
# for set in ko cog; do
#   cp $parsedir/5a_metagenome_key.$set.hierarchy.gexf $heritdir/9l_${set}_hierarchy.gexf
#   cp $parsedir/5c_annotations.$set.enrichment.reformatted.txt $heritdir/9l_${set}_enrichment.aea.txt
#   cp $parsedir/5d_fisher_exact.$set.pvals.txt $heritdir/9l_${set}_enrichment.fisher_exact.txt

# #     # Make a grpah file for tweaking in Gephi for publication
#   python3 $scriptdir/9m_MakeMetagenomeHeritGraph.py -i $heritdir/9k_metagenome_data.herits.txt -g $heritdir/9l_${set}_hierarchy.gexf -o $heritdir/9m_${set}_herit \
#     --pval-cutoff $cutoff --max-herit $max_h2 --fisher $heritdir/9l_${set}_enrichment.fisher_exact.txt
  
#   break
# done

# Make cobined table of category enrichment
python3 $scriptdir/9n_CombineEnrichmentStats.py --aea $heritdir/9l_*_enrichment.aea.txt --fisher $heritdir/9l_*_enrichment.fisher_exact.txt --remove-aea-pval-1 \
  -o $heritdir/9l_metagenome_category_enrichment.combined_pvals.txt

# # # # # # TODO: Break into KEGG / COG terms and graph separately for supplemental data
# # # # # for set in ko cog; do
# # # # #   python3 $scriptdir/9k_PlotMetagenomeHeritData.py --herits $heritdir/9k_metagenome_data.herits.txt  --perms $heritdir/9k_metagenome_data.perms.txt -g $heritdir/9k_metagenome_heritability_pretty.$set \
# # # # #    --pval-cutoff $cutoff -n 115 --figsize 8 10 --key $parsedir/5a_metagenome_key.txt --max-h2 0.98 --traits $heritdir/9j_$set.traits.txt #--debug 
# # # # # #  break
# # # # # done

# TODO: Need to check out the traits with h2 above 0.98 to see if they really are artifacts












####
# Temporary - Collapse things to Genus level for Sahar
# filter_samples_from_otu_table.py -i $qiime_prefix.biom --min_count 10000 -o tmp_sahar/filtered.biom
# summarize_taxa.py -i tmp_sahar/filtered.biom --level 6 --mapping $qiime_prefix.key.tsv -o tmp_sahar