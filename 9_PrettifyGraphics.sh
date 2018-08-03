#! /bin/bash

# This script is not technically part of the analysis pipeline. Instead, it takes the results from the pipeline and formats them into publication-ready graphics, extracts data for tables, etc.

# Software
KRONA=$HOME/Software/Metagenomics/KronaTools-2.7/bin/ktImportText
TASSEL5="perl /home/jgwall/Software/TASSEL/tassel-5-standalone/run_pipeline.pl -Xms10g -Xmx40g"

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
gwasoutdir=$plotdir/9n_GWAS
expressdir=$plotdir/9w_Expression

if [ ! -e $plotdir ]; then mkdir $plotdir; fi
if [ ! -e $kronadir ]; then mkdir $kronadir; fi
if [ ! -e $coredir ]; then mkdir $coredir; fi
if [ ! -e $pairdir ]; then mkdir $pairdir; fi
if [ ! -e $alphadir ]; then mkdir $alphadir; fi
if [ ! -e $betadir ]; then mkdir $betadir; fi
if [ ! -e $heritdir ]; then mkdir $heritdir; fi
if [ ! -e $gwasoutdir ]; then mkdir $gwasoutdir; fi
if [ ! -e $expressdir ]; then mkdir $expressdir; fi

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

# Make alpha diversity plots showing increased diversity with date, subpop, and collector
# targetdir=$alphadir/9j_qiime_diversity_normalized_data/arare_max$rarefaction_level
# make_rarefaction_plots.py -i $targetdir/alpha_div_collated/ -m $qiime_prefix.key.tsv -o $targetdir/alpha_rarefaction_plots --colorby date,subpop,collector --generate_average_tables
# python3 $scriptdir/9j_PlotAlphaDiversity.py -i $targetdir/alpha_rarefaction_plots/average_tables/*.txt -o $alphadir/9j_alpha_diversity.svg

# # Do a statistical test of significance on the above alpha diversity
# targetdir=$alphadir/9j_qiime_diversity_normalized_data/arare_max$rarefaction_level
# pval_cutoff=0.05
# for infile in $targetdir/alpha_div_collated/*.txt; do
#   base=`basename $infile`
#   base=${base/.txt/}
#   outdir=$targetdir/${base}_stats
#   compare_alpha_diversity.py -i $infile -m $qiime_prefix.key.tsv -o $outdir --test_type nonparametric --num_permutations 10000 --categories date,subpop,collector
#   
#   for pvals in $outdir/*_stats.txt; do
#     outfile=${pvals/.txt/.trimmed.txt}
#     Rscript -e "x=read.delim('$pvals'); x=subset(x, x\$p.value <= $pval_cutoff); write.table(x, file='$outfile', sep='\t', quote=F, row.names=F, col.names=T)"
# #     break
#   done
# #   break
# done

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

# # # Look for correlations between beta diversity and different taxa
# Rscript $scriptdir/9i_CorrelateBetaDiversity.r -a $gwasdir/4m_diversity.txt -b $gwasdir/4m_otus.txt -o $betadir/9i_diversity_correlations.txt


# # Plot different PCs on axes and color by different taxonomic groups
# for type in png svg; do
# 
#   # PC1 by methylobacteria
#   python3 $scriptdir/9i_PlotTraitsAndColorize.py -i $gwasdir/4m_diversity.txt $gwasdir/4m_otus.txt -x weighted_unifrac_PC1 -y weighted_unifrac_PC2 \
#     --color qiime.Bacteria.Proteobacteria.Alphaproteobacteria.Rhizobiales.Methylobacteriaceae -o $betadir/9i_weighted_unifrac.pc1.pc2.methylobacteriaceae.$type --cmap Blues
#   
#   # PC2 by bacteroidetes
#   python3 $scriptdir/9i_PlotTraitsAndColorize.py -i $gwasdir/4m_diversity.txt $gwasdir/4m_otus.txt -x weighted_unifrac_PC1 -y weighted_unifrac_PC2 \
#     --color qiime.Bacteria.Bacteroidetes -o $betadir/9i_weighted_unifrac.pc1.pc2.bacteroidetes.$type --cmap Purples
#   
#   # PC3 by sphingomonadaceae
#   python3 $scriptdir/9i_PlotTraitsAndColorize.py -i $gwasdir/4m_diversity.txt $gwasdir/4m_otus.txt -x weighted_unifrac_PC3 -y weighted_unifrac_PC4 \
#     --color qiime.Bacteria.Proteobacteria.Alphaproteobacteria.Sphingomonadales.Sphingomonadaceae -o $betadir/9i_weighted_unifrac.pc3.pc4.sphingomonadaceae.$type --cmap Greens
# 
# # Opted not to plot the following ones
# # # #   # PC1 by rhizobiales
# # # #   python3 $scriptdir/9i_PlotTraitsAndColorize.py -i $gwasdir/4m_diversity.txt $gwasdir/4m_otus.txt -x weighted_unifrac_PC1 -y weighted_unifrac_PC2 \
# # # #     --color qiime.Bacteria.Proteobacteria.Alphaproteobacteria.Rhizobiales -o $betadir/9i_weighted_unifrac.pc1.pc2.rhizobiales.$type
# # # #   # PC3 by gammaproteobacteria
# # # #   python3 $scriptdir/9i_PlotTraitsAndColorize.py -i $gwasdir/4m_diversity.txt $gwasdir/4m_otus.txt -x weighted_unifrac_PC3 -y weighted_unifrac_PC4 \
# # # #     --color qiime.Bacteria.Proteobacteria.Gammaproteobacteria -o $betadir/9i_weighted_unifrac.pc3.pc4.gammaproteobacteria.$type 
# # # # 
# # # #   # PC4 by Hymenobacter
# # # #   python3 $scriptdir/9i_PlotTraitsAndColorize.py -i $gwasdir/4m_diversity.txt $gwasdir/4m_otus.txt -x weighted_unifrac_PC3 -y weighted_unifrac_PC4 \
# # # #     --color qiime.Bacteria.Bacteroidetes.Cytophagia.Cytophagales.Cytophagaceae.Hymenobacter -o $betadir/9i_weighted_unifrac.pc3.pc4.hymenobacter.$type
# # # #    
# # # #    # PC5 by Actinobacteria
# # # #   python3 $scriptdir/9i_PlotTraitsAndColorize.py -i $gwasdir/4m_diversity.txt $gwasdir/4m_otus.txt -x weighted_unifrac_PC4 -y weighted_unifrac_PC5 \
# # # #     --color qiime.Bacteria.Actinobacteria -o $betadir/9i_weighted_unifrac.pc4.pc5.actinobacteria.$type
#   
# #   break
# done

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
# metadir=$heritdir/9j_metagenome_tables
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

# # Plot heritabilities 
# cutoff=0.001
# # OTUs
# python3 $scriptdir/9k_PlotHeritability_publication.py -i $gwasdir/4s_narrow_herit_bigperms_otus/*/4r_*_stats.txt -g $heritdir/9k_otu_heritability_pretty \
#   --pval-cutoff $cutoff -n 20 --figsize 5.5 6 #--debug

# Separate heritability plotting into 2 steps so don't have to load everything each time
# python3 $scriptdir/5_AssembleMetagenomeHeritData.py -i $gwasdir/4s_narrow_herit_bigperms_metagenome*/*/4r_*_stats.txt -o $heritdir/9k_metagenome_data #--debug


###########
# Metagenome heritability graphs, etc.
###########

# # Plot heritability of  meteagenome terms
# cutoff=0.001
# for set in ko cog; do
#   python3 $scriptdir/9k_PlotMetagenomeHeritData.py --herits $heritdir/9k_metagenome_data.herits.txt  --perms $heritdir/9k_metagenome_data.perms.txt -g $heritdir/9k_metagenome_heritability_pretty.$set \
#    --pval-cutoff $cutoff -n 115 --figsize 8 10 --key $parsedir/5a_metagenome_key.txt --max-h2 0.98 --traits $heritdir/9j_$set.traits.txt #--debug 
# #  break
# done

# # Plot heritability of top meteagenome terms, flagging ones that are likely artifacts
# cutoff=0.001
# blups=$gwasdir/4t_combined_good_traits.tassel.txt
# python3 $scriptdir/9k_PlotArtifactualMetagenomeHeritData.py --herits $heritdir/9k_metagenome_data.herits.txt  --perms $heritdir/9k_metagenome_data.perms.txt -g $heritdir/9k_metagenome_heritability.artifacts \
#   --pval-cutoff $cutoff --num-goodtraits 8 --figsize 10 6 --key $parsedir/5a_metagenome_key.txt --max-h2 0.98 --seed 2 --blups $blups #--debug 



# Make metagenome heritability graph
# cutoff=0.001
# max_h2=0.98 # Things above this are thought to be artifacts due to their distribution of permuted heritabilities
# for set in cog ko; do
#   cp $parsedir/5a_metagenome_key.$set.hierarchy.gexf $heritdir/9l_${set}_hierarchy.gexf
#   cp $parsedir/5c_annotations.$set.enrichment.reformatted.txt $heritdir/9l_${set}_enrichment.aea.txt
#   cp $parsedir/5d_fisher_exact.$set.pvals.txt $heritdir/9l_${set}_enrichment.fisher_exact.txt
# 
#     # Make a grpah file for tweaking in Gephi for publication
#   python3 $scriptdir/9m_MakeMetagenomeHeritGraph.py -i $heritdir/9k_metagenome_data.herits.txt -g $heritdir/9l_${set}_hierarchy.gexf -o $heritdir/9m_${set}_herit \
#     --pval-cutoff $cutoff --max-herit $max_h2 --fisher $heritdir/9l_${set}_enrichment.fisher_exact.txt
#   
# #   break
# done

# # Make cobined table of category enrichment
# python3 $scriptdir/9n_CombineEnrichmentStats.py --aea $heritdir/9l_*_enrichment.aea.txt --fisher $heritdir/9l_*_enrichment.fisher_exact.txt --remove-aea-pval-1 \
#   -o $heritdir/9l_metagenome_category_enrichment.combined_pvals.txt


############
# GWAS Analysis
############

# # Get a list of all significant hits for all traits; used to filter down to just the right hit-trait combinations (but want to work with a much smaller file than the raw data)
# ls $gwasdir/4u_gwas/4u_orig_gwas | sed -r "s/_chr.+//" | uniq | sort | uniq > $gwasoutdir/9n_hits.all_traits.txt
# Rscript -e "x=read.delim('$gwasdir/4x_gwas_hits.clustered.best.txt', stringsAsFactors=F); write(unique(x\$Marker), file='$gwasoutdir/9n_hits.markers.txt');" \
#   -e "write(unique(x\$Trait), file='$gwasoutdir/9n_hits.traits.txt')"
# zcat $gwasdir/4u_gwas/4u_orig_gwas/*sitefile.txt.gz | head -n 1 >  $gwasoutdir/9n_hits.sitefiles.txt 
# find $gwasdir/4u_gwas/4u_orig_gwas/ -name "*sitefile.txt.gz" | parallel zcat | grep --file $gwasoutdir/9n_hits.markers.txt >> $gwasoutdir/9n_hits.sitefiles.txt 

# # # Get a list of "bad" traits whose heritabilities are too high and are likely artifacts and filter them out from the clustered trait analysis
# max_h2=0.98
# max_pval=0.001
# python3 $scriptdir/9n_FilterClusteredHits.py -i $gwasdir/4x_gwas_hits.clustered.best.txt --herits $heritdir/9k_metagenome_data.herits.txt --max-herit $max_h2 \
#   -o $gwasoutdir/9n_gwas_hits.clustered.best.goodtraits.txt -t $gwasoutdir/9n_hits.badtraits.txt
# Rscript -e "all=scan('$gwasoutdir/9n_hits.all_traits.txt', what=character()); bad=scan('$gwasoutdir/9n_hits.badtraits.txt', what=character())" \
#   -e "excluded=intersect(all, bad); write(excluded, file='$gwasoutdir/9n_hits.excluded_traits.txt')"

# # Same as above, but using all the hits instead of just the best one per cluster; this is what was used for a Supplemental Table
# python3 $scriptdir/9n_FilterClusteredHits.py -i $gwasdir/4x_gwas_hits.clustered.all.txt --herits $heritdir/9k_metagenome_data.herits.txt --max-herit $max_h2 \
#   -o $gwasoutdir/9n_gwas_hits.clustered.all.goodtraits.txt -t $gwasoutdir/9n_hits.all.badtraits.txt
# Rscript -e "all=scan('$gwasoutdir/9n_hits.all_traits.txt', what=character()); bad=scan('$gwasoutdir/9n_hits.all.badtraits.txt', what=character())" \
#   -e "excluded=intersect(all, bad); write(excluded, file='$gwasoutdir/9n_hits.all.excluded_traits.txt')"


# # Identify traits with the most number of hits below a given p-value
# for p in 0.05 0.01; do
#   python3 $scriptdir/9o_TabulateGwasHits.py -i $gwasoutdir/9n_gwas_hits.clustered.best.goodtraits.txt -p $p -o $gwasoutdir/9o_hits_per_trait.p$p.txt
# #   python3 $scriptdir/9o_GetGwasSitefileData.py -i $gwasoutdir/9n_gwas_hits.clustered.best.goodtraits.txt -p $p --sitefile $gwasoutdir/9n_hits.sitefiles.txt  -o $gwasoutdir/9o_sitefile_data.p$p.txt
# #   break
# done

# Plot Manhatten plots of flowering time + 2 others (unifract not included in publication; used for a seminar)
# for trait in flowering_time log_COG0181_1489 log_K02769_7737; do
#   python3 $scriptdir/9p_GatherGwasData.py --raw $gwasdir/4u_gwas/4u_orig_gwas/${trait}_chr*.pvals.txt --perms $gwasdir/4u_gwas/4v_gwas_results/${trait}.pvals.txt -o $gwasoutdir/9p_$trait.pvals.txt
# #   break
# done
# python3 $scriptdir/9p_PlotGwasData.py -i $gwasoutdir/9p_*.pvals.txt --chromlengths $scriptdir/0_Zmays_agpv3_chrom_lengths.txt -o $gwasoutdir/9p_manhatten_plots \
#   --flowering-genes $scriptdir/0_flowering_candidate_genes.csv --candidate-window 1000000  #--sparsify 0.01

# # # Get distances between each hit and the closest flowering time candidate
# Rscript -e "x=read.delim('$gwasoutdir/9n_gwas_hits.clustered.best.goodtraits.txt'); x=subset(x, x\$Trait=='flowering_time')"\
#   -e "write.table(x, file='$gwasoutdir/9q_flowering_time.clusters.txt', sep='\t', quote=F, row.names=F, col.names=T)"
# python3 $scriptdir/9q_GetNearestFloweringCandidates.py -i $gwasoutdir/9p_flowering_time.pvals.txt --flowering-genes $scriptdir/0_flowering_candidate_genes.csv \
#   -o $gwasoutdir/9q_flowering_time.all.nearest_candidates.txt    # All hits
# python3 $scriptdir/9q_GetNearestFloweringCandidates.py -i $gwasoutdir/9q_flowering_time.clusters.txt --flowering-genes $scriptdir/0_flowering_candidate_genes.csv \
#   -o $gwasoutdir/9q_flowering_time.clusters.nearest_candidates.txt  # just the best hit in each cluster


# # Get the closest genes for each hit. Cutoff of 0.01 is for text b/c those are the believable ones (low pval and best in cluster); cutoff of 0.1 is for supplemental so people can see everything.
# genes=$scriptdir/0_agpv3_genes.gff
# python3 $scriptdir/9q_GetNearestCandidateGenes.py -i $gwasoutdir/9n_gwas_hits.clustered.best.goodtraits.txt -o $gwasoutdir/9q_all_traits.nearest_candidates --pval-cutoff 0.01 --gff $genes \
#   --key $heritdir/9j_metagenome_key.combined.txt #--debug
# python3 $scriptdir/9q_GetNearestCandidateGenes.py -i $gwasoutdir/9n_gwas_hits.clustered.all.goodtraits.txt -o $gwasoutdir/9q_all_traits.nearest_candidates.pval_0.1 --pval-cutoff 0.1 --gff $genes \
#   --key $heritdir/9j_metagenome_key.combined.txt #--debug


# # # Table of regions that have repeated hits
# metagenome_key=$parsedir/5a_metagenome_key.txt
# for winsize in 1000000 100000 10000; do
#   python3 $scriptdir/9r_TabulateHitClusters.py -i $gwasoutdir/9n_gwas_hits.clustered.best.goodtraits.txt -o $gwasoutdir/9r_common_hits.winsize$winsize --winsize $winsize \
#     --cors $parsedir/5f_blup_correlations.txt --min-score 1.5 --pval-cutoff 0.01
#   python3 $scriptdir/9s_PrettifyShortHitTable.py -i $gwasoutdir/9r_common_hits.winsize$winsize.windows_short.filtered.txt -m $metagenome_key -o $gwasoutdir/9s_clusters_short.win$winsize.pretty.txt
# #   break
# done

# # Just anlayze 100k in depth
# winsize=100000
# blups=$gwasdir/4u_combined_good_traits.update_flowering.tassel.txt
# # Plot correlations among traits
# python3 $scriptdir/9t_AnalyzeClusteredTraits.py -i $gwasoutdir/9r_common_hits.winsize$winsize.windows_short.filtered.txt -o $gwasoutdir/9t_common_hits.winsize$winsize.trait_correlations --blups $blups

# # Look at linkage disequilibrium in this region
# genos=$gwasdir/4u_gwas/all_genos.hmp.txt.gz
# zcat $genos | tail -n +2 | cut -f1-4 > $gwasoutdir/9u_ld.all_markers.txt
# python3 $scriptdir/9u_GetSitesForLd.py -i $gwasoutdir/9s_clusters_short.win$winsize.pretty.txt -m $gwasoutdir/9u_ld.all_markers.txt -o $gwasoutdir/9u_ld.target_markers.txt
# grep -e "gene" -e "exon" $HOME/Projects/0_RawData/MaizeGenome/Zea_mays.AGPv3.30.gff3 > $gwasoutdir/9u_maize_genes.gff3

# while read marker sitenum; do
#   $TASSEL5 -h $genos -ld -ldType SiteByAll -ldTestSite $sitenum -ldHetTreatment Homozygous -export $gwasoutdir/9u_ld.$marker.txt
#   # Plot at two different scales
#   python3 $scriptdir/9u_PlotLd.py -i $gwasoutdir/9u_ld.$marker.txt -o $gwasoutdir/9u_ld.$marker.ld_window --big-winsize 10000000 --small-winsize 100000 --gff $gwasoutdir/9u_maize_genes.gff3
# #   break
# done < $gwasoutdir/9u_ld.target_markers.txt
# 
# 
# # # Now look at phenotype distributions by alleles to make sure that separation is real and not some strange artifact
# genos=$gwasdir/4u_gwas/all_genos.hmp.txt.gz
# blups=$gwasdir/4u_combined_good_traits.update_flowering.tassel.txt
# $TASSEL5 -h $genos -includeSiteNamesInFile $gwasoutdir/9n_hits.markers.txt -export $gwasoutdir/9v_hit_genos.hmp.txt
# python3 $scriptdir/9v_PlotPhenosByAllelesForHits.py -i $gwasoutdir/9n_gwas_hits.clustered.best.goodtraits.txt --genos $gwasoutdir/9v_hit_genos.hmp.txt --phenos $blups --pval-cutoff 0.01 -o $gwasoutdir/9v_hit_distributions.png




##############
# Expression correlation
##############

# python3 $scriptdir/9w_CorrelateTraitsWithExpression.py -i $gwasdir/4m_diversity.txt --traits weighted_unifrac_PC1 --genes GRMZM2G031545 --expression $scriptdir/0_kremling_expression_data.txt \
#   -o $expressdir/9w_diversity.expression.png --key $scriptdir/0_kremling_expression_key.txt
# python3 $scriptdir/9w_CorrelateTraitsWithExpression.py -i $gwasdir/4m_metagenome.txt --traits D_Arginine_and_D_ornithine_metabolism_11861 K14977_11335 COG2319_3037 COG3484_4168 K07395_9407 --genes GRMZM2G031545 \
#   --expression $scriptdir/0_kremling_expression_data.txt -o $expressdir/9w_metagenome.expression.png --key $scriptdir/0_kremling_expression_key.txt
# 
# # SVG for supplemental figure
# python3 $scriptdir/9w_CorrelateTraitsWithExpression.py -i $gwasdir/4m_metagenome.txt --traits D_Arginine_and_D_ornithine_metabolism_11861 K14977_11335 COG2319_3037 COG3484_4168 K07395_9407 --genes GRMZM2G031545 \
#   --expression $scriptdir/0_kremling_expression_data.txt -o $expressdir/9w_metagenome.expression.svg --key $scriptdir/0_kremling_expression_key.txt 

# # Quick check that no correlation exists with other nearby genes; first assemble a list, then pass it all to the plot script
# cut -f9 $gwasoutdir/9u_ld.S7_165387694.ld_window.genes.gff | grep 'ID=gene:' | sed -r "s|ID=gene:([^;]+);.+|\1|"> $expressdir/9x_chr7_target_genes.txt
# readarray genes < $expressdir/9x_chr7_target_genes.txt
# genes=`echo ${genes[@]}`
# python3 $scriptdir/9w_CorrelateTraitsWithExpression.py -i $gwasdir/4m_metagenome.txt --traits D_Arginine_and_D_ornithine_metabolism_11861 K14977_11335 COG2319_3037 COG3484_4168 K07395_9407 --genes $genes \
#   --expression $scriptdir/0_kremling_expression_data.txt -o $expressdir/9x_metagenome.expression.png --key $scriptdir/0_kremling_expression_key.txt 


###############
# Minor cleanup stuff
###############

# # This isn't a graphic per se, just a place to collect all the packages and libraries I need to cite in the paper
# grep --no-filename "import" $scriptdir/*.py | sed -e "s/import //" -e "s/from //" | sort | uniq  > $plotdir/9z_python_imports.txt
# grep --no-filename "library(" $scriptdir/*.r | sort | uniq > $plotdir/9z_r_libraries.txt

# # Collect file names for the SRA submission
# Rscript $scriptdir/9z_CollectFileNamesForSra.r -i $rawdir/*.fastq.gz -o $plotdir/9z_sra_file_names.txt