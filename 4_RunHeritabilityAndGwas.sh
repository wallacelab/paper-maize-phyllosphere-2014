#! /bin/bash 

# Collapse OTU counts to BLUPs and then run statistics off of those

TASSEL5="perl $HOME/Software/TASSEL/tassel-5-standalone/run_pipeline.pl -Xms10g -Xmx40g"

scriptdir=$1
gwasdir=$2
qiime_prefix=$3
rarefaction_level=$4
min_counts_across_samples=$5
num_pcs=$6
flowering_data=$7
cofactors="$8"
genos=$9
aliasfile=${10}
top_n=${11}
nperms=${12}
bigperms=${13}
perm_cutoff=${14}
gwas_perms=${15}
chromlengths=${16}
TASSEL5="${17}"
cogdir=${18}
ncores=${19}
gwas_procs=${20}
empirical_pval_cutoff=${21}
all_flowering=${22}
cluster_cutoff=${23}

biom=$qiime_prefix.biom
keyfile=$qiime_prefix.key.tsv
taxonomy=$qiime_prefix.taxonomy.txt

###############
# Get basic stats on biom file
###############

# # Get a summary of the top N OTUs in each. (Useful for deciding what value of top_n to pass)
# otu_summary=$gwasdir/4a_otu_counts.txt
# biom summarize-table -i $biom -o $otu_summary --observations
# Rscript $scriptdir/4a_BiomStats.r -i $otu_summary -o $gwasdir/4a_otu_cutoffs.txt --outgraph $gwasdir/4a_otu_cutoffs.png --top-n $top_n --topfile $gwasdir/4a_top${top_n}_otus.txt

# # Do the same for samples
# samplesummary=$gwasdir/4b_sample_counts.txt
# samplecut=$gwasdir/4b_sample_cutoffs.txt
# samplegraph=$gwasdir/4b_sample_cutoffs.png
# biom summarize-table -i $biom -o $samplesummary 
# Rscript $scriptdir/4a_BiomStats.r -i $samplesummary -o $samplecut --outgraph $samplegraph --log xy


###############
# Prepare OTUs and higher-level taxa
###############

# #   Collapse individual OTUs into higher-level taxonomy
# biom convert -i $biom -o $gwasdir/4c_otus.txt --to-tsv --table-type "OTU table"
# python3 $scriptdir/4d_CollapseOtusByTaxonomy.py -i $gwasdir/4c_otus.txt -o $gwasdir/4d_collapsed_taxa.txt -t $taxonomy --prefix qiime

# # Normalize samples to [rarefaction_level] reads (based on original counts), remove any samples with less than that, and remove any OTUs/taxa that have counts in fewer than [min_counts_across_samples] samples
# Rscript $scriptdir/4e_NormalizeAndFilterTaxa.r -i $gwasdir/4d_collapsed_taxa.txt -o $gwasdir/4e_collapsed_taxa.filtered.txt -p $gwasdir/4e_collapsed_taxa.filtered.png \
#    --normalize-count $rarefaction_level --min-count-across-samples $min_counts_across_samples --sample-summary $samplesummary


################
# Prepare diversity metrics
################

# # Normalize data to [rarefaction_level]; essentially the asymptotic result of an infinite number of rarefactions. Min-count set low enough only things with all zeroes should be removed
# Rscript $scriptdir/4e_NormalizeAndFilterTaxa.r -i $qiime_prefix.txt -o $gwasdir/4f_orig_otus.filtered.txt -p $gwasdir/4f_orig_otus.filtered.png \
#    --normalize-count $rarefaction_level --min-count-across-samples 0.0000001 --sample-summary $samplesummary --skip 1  

# # Convert to a real biom file and add taxonomy back in
# biom convert -i $gwasdir/4f_orig_otus.filtered.txt -o $gwasdir/4f_orig_otus.filtered.biom --to-hdf5
# biom add-metadata -i  $gwasdir/4f_orig_otus.filtered.biom -o $gwasdir/4f_orig_otus.filtered.with_metadata.biom --observation-metadata-fp $qiime_prefix.taxonomy.txt \
#   --sc-separated taxonomy --observation-header OTUID,taxonomy


# # Alpha diversity metrics
# metrics=$(alpha_diversity.py --show_metrics | grep "^Known" | sed -e "s/Known metrics are: //" -e "s/ //g" -e "s/michaelis_menten_fit,//" )  # Michaelis-Menten hangs analysis, so skip
# alpha_diversity.py -i $gwasdir/4f_orig_otus.filtered.with_metadata.biom -o $gwasdir/4g_orig_otus.filtered.alpha_diversity.txt --metrics $metrics --tree_path $qiime_prefix.tre

# # Beta diversity metrics
# metrics="unweighted_unifrac,weighted_unifrac,bray_curtis"
# betadir=$gwasdir/4g_beta_diversity
# pcdir=$gwasdir/4g_principal_coordinates
# beta_diversity.py -i $gwasdir/4f_orig_otus.filtered.with_metadata.biom -o $betadir --metrics $metrics --tree_path $qiime_prefix.tre
# principal_coordinates.py -i $betadir -o $pcdir
# python3 $scriptdir/4g_ConvertQiimePcsToPlainMatrix.py -i $pcdir/*.txt -o $gwasdir/4g_orig_otus.filtered.beta_diversity.txt -g $gwasdir/4g_orig_otus.filtered.beta_diversity.histograms.png \
#     --num-pcs $num_pcs --screeplot $gwasdir/4g_orig_otus.filtered.beta_diversity.scree.png

# # Merge alpha and beta diversity into a single file and format similar to otus and metagenome (= samples in columns and metrics in rows)
# Rscript -e "alpha=read.delim('$gwasdir/4g_orig_otus.filtered.alpha_diversity.txt', row.names=1); beta=read.delim('$gwasdir/4g_orig_otus.filtered.beta_diversity.txt', row.names=1)" \
#   -e "output=merge(alpha, beta, by='row.names'); write.table(t(output), file='$gwasdir/4g_orig_otus.filtered.both_diversity.txt', sep='\t', quote=F, col.names=F)"


################
# Prepare imputed metagenome from same rarefied-and-filtered set as above
################

# # Trim rarefied QIIME OTU table down to just reference OTUs.
# grep -e "New.Reference" -e "New.CleanUp.Reference" $qiime_prefix.observations.txt | sed -r "s/:.+//" > $gwasdir/4h_qiime.nonreference_otus.txt
# filter_otus_from_otu_table.py -i $gwasdir/4f_orig_otus.filtered.with_metadata.biom -o $gwasdir/4h_closed_reference_otus.biom --otu_ids_to_exclude_fp $gwasdir/4h_qiime.nonreference_otus.txt
# biom summarize-table -i $gwasdir/4h_closed_reference_otus.biom -o $gwasdir/4h_closed_reference_otus.samples.txt  
# biom summarize-table -i $gwasdir/4h_closed_reference_otus.biom -o $gwasdir/4h_closed_reference_otus.observations.txt  --observations

# # Use PICRUST to infer metagenome content 
# normalize_by_copy_number.py -i $gwasdir/4h_closed_reference_otus.biom -o $gwasdir/4i_closed_reference_otus.fix_copy_number.biom
# for set in ko cog; do
#   echo "Predicting metagenome on QIIME samples for prediction set '$set'"
#   predict_metagenomes.py -i $gwasdir/4i_closed_reference_otus.fix_copy_number.biom -o $gwasdir/4j_predicted_metagenome.$set.biom \
#     -a $gwasdir/4j_predicted_metagenome.$set.nsti_vals.txt --type_of_prediction $set 
# done
# 
# # Collapse individual KO and COG terms by levels. (KEGG has 3, COG has 2)
# for level in 1 2 3; do
#   categorize_by_function.py -i $gwasdir/4j_predicted_metagenome.ko.biom -c KEGG_Pathways -l $level -o $gwasdir/4j_predicted_metagenome.ko.l$level.biom
# done
# for level in 1 2; do
#   categorize_by_function.py -i $gwasdir/4j_predicted_metagenome.cog.biom -c COG_Category -l $level -o $gwasdir/4j_predicted_metagenome.cog.l$level.biom
# done
# 
# # Convert to text
# for infile in $gwasdir/4j_*.biom; do
#   outfile=${infile/.biom/.biom.txt}
#   biom convert -i $infile -o $outfile --to-tsv
# done

# # Merge all metagenome tables into one
# python3 $scriptdir/4k_CombineOtuTables.py -i $gwasdir/4j_*.biom.txt -o $gwasdir/4k_combined_metagenome.txt --skip 1 --make-unique

# # Remove predicted functions that aren't present in enough samples; no actual normalization done, just filtering
# Rscript $scriptdir/4e_NormalizeAndFilterTaxa.r -i $gwasdir/4k_combined_metagenome.txt -o $gwasdir/4l_predicted_metagenome.filtered.txt -p $gwasdir/4l_predicted_metagenome.filtered.png \
#   --min-count-across-samples $min_counts_across_samples --sample-summary $samplesummary --skip-normalization  --normalize-count 1



######################
# Get BLUPs for each dataset
######################

# # Copy everything over to get into the same format
# python3 $scriptdir/4m_StandardizeDataTables.py -i $gwasdir/4e_collapsed_taxa.filtered.txt -o $gwasdir/4m_otus.txt --fix-traitnames
# python3 $scriptdir/4m_StandardizeDataTables.py -i $gwasdir/4g_orig_otus.filtered.both_diversity.txt -o $gwasdir/4m_diversity.txt --fix-traitnames
# python3 $scriptdir/4m_StandardizeDataTables.py -i $gwasdir/4l_predicted_metagenome.filtered.txt -o $gwasdir/4m_metagenome.txt --fix-traitnames

# # Make a file of just flowering time
# python3 $scriptdir/4m_StandardizeDataTables.py -i $gwasdir/4m_otus.txt -o $gwasdir/4m_flowering.txt -f $flowering_data --jitter --skip-orig-data

# # Identify significant cofactors; done first time to determine what to include in next step
# for set in diversity otus metagenome; do
#   Rscript $scriptdir/4n_FindSignificantCofactors.r -i $gwasdir/4m_$set.txt -k $keyfile -o $gwasdir/4n_cofactor_pvals.$set.1.txt  \
#     --cofactors-to-test tissue_date dna_plate collector --num-cores $ncores --plotfile $gwasdir/4n_cofactor_pvals.$set.1.png #--debug
#   Rscript $scriptdir/4n_FindSignificantCofactors.r -i $gwasdir/4m_$set.txt -k $keyfile -o $gwasdir/4n_cofactor_pvals.$set.2.txt  \
#     --covariates collector --cofactors-to-test tissue_date dna_plate --num-cores $ncores --plotfile $gwasdir/4n_cofactor_pvals.$set.2.png #--debug
#   Rscript $scriptdir/4n_FindSignificantCofactors.r -i $gwasdir/4m_$set.txt -k $keyfile -o $gwasdir/4n_cofactor_pvals.$set.3.txt  \
#     --covariates tissue_date --cofactors-to-test collector dna_plate --num-cores $ncores --plotfile $gwasdir/4n_cofactor_pvals.$set.3.png #--debug
#   Rscript $scriptdir/4n_FindSignificantCofactors.r -i $gwasdir/4m_$set.txt -k $keyfile -o $gwasdir/4n_cofactor_pvals.$set.4.txt  \
#     --covariates collector tissue_date --cofactors-to-test dna_plate --num-cores $ncores --plotfile $gwasdir/4n_cofactor_pvals.$set.4.png  #--debug
# #   break
# done

# # Split metagenome into multiple susbets because it's so large
# splits=8
# Rscript -e "x=read.delim('$gwasdir/4m_metagenome.txt'); x=split(x, rep(1:$splits, nrow(x))[1:nrow(x)]);" \
#   -e "lapply(1:$splits, function(i){write.table(x[[i]], file=paste('$gwasdir/4m_metagenome',i,'.txt', sep=''), sep='\t', quote=F, row.names=F, col.names=T)})"

# # Make BLUPs to combine day/night and biological replicates into single samples per genotype
# for set in diversity otus metagenome1 metagenome2 metagenome3 metagenome4 metagenome5 metagenome6 metagenome7 metagenome8; do
#   Rscript $scriptdir/4o_CalculateBlups.r -i $gwasdir/4m_$set.txt -k $keyfile -o $gwasdir/4o_$set.blups.tassel.txt \
#     --target-column Description --covariates $cofactors --num-cores $ncores --plotfile $gwasdir/4o_$set.blups.distributions.png --num-traits-to-plot 40 --log-compare # --debug
# #   break
# done

# # Do flowering separately with a different model (since it shouldn't be affected by the covariates, although it is partially confounded with them)
# Rscript $scriptdir/4o_CalculateBlups.r -i $gwasdir/4m_flowering.txt -k $keyfile -o $gwasdir/4o_flowering.blups.tassel.txt \
#     --target-column Description --num-cores $ncores --plotfile $gwasdir/4o_flowering.blups.distributions.png --num-traits-to-plot 1 --log-compare # --debug

##################
# Narrow-sense heritability
##################

# # Match up sample names so match the ones in the genotype file
# # NOTE: 6 CML lines are not present in the genotype file (CML411, CML504, CML505, CML84, CML85, CML96) and do not appear in the public Panzea GBS dataset
# python3 $scriptdir/4p_MatchGenoNamesToSamples.py -m $keyfile -g  $genos -o $gwasdir/4p_samples_genos_key.txt -a $aliasfile
# for infile in $gwasdir/4o_*.blups.tassel.txt; do
#   python3 $scriptdir/4p_CorrectGenoNamesInBlups.py -i $infile -o ${infile/4o_/4p_} -k $gwasdir/4p_samples_genos_key.txt 
# #   break
# done

# # Files for narrow-sense heritability
# kinship=${genos/.hmp.txt.gz/.kinship.txt}
# dummy_genos=$gwasdir/4q_dummy_genos.hmp.txt
# flowering_covariate=$gwasdir/4p_flowering_covariate.tassel.txt
# flowering_blups=$gwasdir/4p_flowering.blups.tassel.txt
# zcat $genos | head -n 2 > $dummy_genos

# # Turn flowering time into a TASSEL covariate file (involves changing just one part of the header)
# sed "s/data$/covariate/" $gwasdir/4p_flowering.blups.tassel.txt > $flowering_covariate

# # Loop over datasets and do narrow-sense heritability
# for blups in $gwasdir/4p_*.blups.tassel.txt; do
#   set=${blups/*4p_/}
#   set=${set/.blups.tassel.txt/}
#   blups_with_flowering=$gwasdir/4q_blups_with_flowering.$set.tassel.txt
# 
#   if [[ $set == 'flowering' ]] ; then continue; fi # No need to do flowering on its own; can fold into everything else.
#   # Add flowering to file
#   $TASSEL5 -fork1 -t $blups -fork2 -t $flowering_blups -combine3 -input1 -input2 -union -export $blups_with_flowering
# 
#   # Directory to store intermediate files (b/c there's a lot of them)
#   narrowdir=$gwasdir/4q_narrow_herit_$set
#   if [ ! -e $narrowdir ]; then mkdir $narrowdir; fi
# 
#   cutoff=0.02
#   outprefix=$narrowdir
#   $scriptdir/4q_RunTasselHeritabilityPerms_WithCovariates.sh $scriptdir $blups_with_flowering $dummy_genos $kinship $nperms $cutoff $narrowdir $outprefix $flowering_covariate flowering_time "$TASSEL5"
# #   break
# done


# # Filter BLUPs for traits that passed initial filtering, then rerun with higher number of permutations
# tail -q -n +2  $gwasdir/4q_*.good_traits.txt | cut  -f1 | sort | uniq > $gwasdir/4r_combined_good_traits.txt
# for blups in $gwasdir/4q_blups_with_flowering.*.tassel.txt; do
#   set=${blups/*4q_blups_with_flowering./}
#   set=${set/.tassel.txt/}
#   filtered_blups=$gwasdir/4r_good_traits.$set.tassel.txt
#   
#   python3 $scriptdir/4r_FilterGoodTraits.py -i $blups --traitfile $gwasdir/4r_combined_good_traits.txt -o $filtered_blups
# 
#   narrowdir=$gwasdir/4s_narrow_herit_bigperms_$set
#   if [ ! -e $narrowdir ]; then mkdir $narrowdir; fi
#   outprefix=$narrowdir
#   $scriptdir/4q_RunTasselHeritabilityPerms_WithCovariates.sh $scriptdir $filtered_blups $dummy_genos $kinship $bigperms $perm_cutoff $narrowdir $outprefix $flowering_covariate flowering_time "$TASSEL5"
  # 
# done 

# # Get KEGG annotations for metagenome
# tail -q -n+2 $gwasdir/4s_narrow_herit_bigperms_metagenome*.good_traits.txt | cut -f1 | sed -r "s/^log_//" | grep "^K" | sed -r "s/_.+//"> $gwasdir/4s_kegg_accessions.txt
# Rscript $scriptdir/4t_RetrieveKeggAnnotations.r -i $gwasdir/4s_kegg_accessions.txt -o $gwasdir/4s_kegg_annotations.txt
# 
# # Collate and parse the metagenome results better; accesses annotations for COG and KEGG data
# python3 $scriptdir/4s_CollateAndParseMetagenomeHerits.py -i $gwasdir/4s_narrow_herit_bigperms_metagenome*.good_traits.txt -o $gwasdir/4s_narrow_herit_bigperms.metagenome_combined.txt --cog-names $cogdir/cognames2003-2014.txt \
#   --cog-functions $cogdir/fun2003-2014.txt --kegg $gwasdir/4s_kegg_annotations.txt
  
# TODO: Identify common pathways in the above?

######
# Run GWAS for the traits that passed filtering, using chromosome-specific kinship matrices to get chromosome-specific residuals, then running these against permuted genotypes to get p-values
######



# # # Cut down and combine BLUP files
# tail -q -n +2  $gwasdir/4s_*.good_traits.txt | cut  -f1 | sort | uniq > $gwasdir/4t_combined_good_traits.txt
# forks=""; runs=""; inputs=""    # Components of TASSEL command to combine phenotypes
# for blups in $gwasdir/4q_blups_with_flowering.*.tassel.txt; do
#   set=${blups/*4q_blups_with_flowering./}
#   set=${set/.tassel.txt/}
#   filtered_blups=$gwasdir/4t_best_traits.$set.tassel.txt
# #   python3 $scriptdir/4r_FilterGoodTraits.py -i $blups --traitfile $gwasdir/4t_combined_good_traits.txt -o $filtered_blups
#   
#   # TASSEL command buildup
#   forks="$forks -fork$set -t $filtered_blups"
#   if [ $set != "otus" ] ; then forks="$forks -excludeLastTrait" ; fi    # Only let 1 flowering time trait come through
#   inputs="$inputs -input$set"
#   runs="$runs -runfork$set"
# done 
# $TASSEL5 $forks -combine000 $inputs -intersect -export $gwasdir/4t_combined_good_traits.tassel.txt $runs    


# # # GWAS parameters
# outdir=$gwasdir/4u_gwas
# permdir=$outdir/4u_permuted_genos
# if [ ! -e $outdir ] ; then mkdir $outdir; fi
# if [ ! -e $permdir ] ; then mkdir $permdir; fi
# 
# gwas_blups=$gwasdir/4t_combined_good_traits.tassel.txt
# kinship=${genos/.hmp.txt.gz/.kinship.txt}


# # Filter genotypes by those in my sample file (to save space), set all hets to missing, split by chromosome, and get kinship matrices including all _but_ that chromosome
# # Note: hets removed because caused issues in downstream mapping due to rare outlier phenotypes with het genos. Minimally affect kinship matrices (>99.9% correlated), so heritability is unaffected.
# cut -f1 $gwas_blups | tail -n+4 > $gwasdir/4t_combined_good_traits.taxa.txt
# $TASSEL5 -h $genos -homozygous -includeTaxaInFile $gwasdir/4t_combined_good_traits.taxa.txt -export $outdir/all_genos.hmp.txt.gz
# for chr in $(seq 1 10); do
#   $TASSEL5 -h $genos -homozygous -includeTaxaInFile $gwasdir/4t_combined_good_traits.taxa.txt -separate $chr -export $outdir/chr$chr.hmp.txt.gz
#   zcat $outdir/chr$chr.hmp.txt.gz | head -n2  | gzip > $outdir/chr$chr.dummy_genos.hmp.txt.gz
# done

# # Make kinship matrices from all other chromosomes
# for chr in $(seq 1 10); do
#   # Build TASSEL arguments
#   forks=""
#   inputs=""
#   runs=""
#   for fork in $(seq 1 10); do
#     if [ $fork -eq $chr ]; then continue; fi # Skip the chromosome in question
#     forks="$forks -fork$fork -h $outdir/chr$fork.hmp.txt.gz"
#     inputs="$inputs -input$fork"
#     runs="$runs -runfork$fork"
#   done
#   #Run
#   $TASSEL5 $forks -combine000 $inputs -union -ck -export $outdir/chr$chr.other_chroms.kinship.txt $runs
# #   break
# done
# 
# Permute genotypes
# for chr in $(seq 1 10); do
#   genos=$outdir/chr$chr.hmp.txt.gz
#   perm_prefix=$permdir/chr$chr.perm
#   python3 $scriptdir/4u_PermuteHapmapGenotypes.py -i $genos -o $perm_prefix --outpostfix .hmp.txt -n $gwas_perms --seed $chr #--debug
#   pigz -f $perm_prefix*.hmp.txt # Parallel zipping much faster that using python internal gzip
# #   break
# done

# # Interpolate flowering covariate so I don't lose 10% of samples due to missing flowering time in this particular field
# Rscript $scriptdir/4u_InterpolateMissingFloweringTime.r -i $gwasdir/4p_flowering_covariate.tassel.txt -f $all_flowering -t $gwasdir/4t_combined_good_traits.taxa.txt \
#   -o $gwasdir/4u_flowering_covariate.interpolated.tassel.txt

# # Replace flowering time in BLUPs with this interpolated one
# Rscript $scriptdir/4u_ReplaceFloweringBlupsWithInterpolated.r -b $gwas_blups -i $gwasdir/4u_flowering_covariate.interpolated.tassel.txt -o $gwasdir/4u_combined_good_traits.update_flowering.tassel.txt

# # Run GWAS (Note: this is actually a fairly involved script in itself, and was run on UGA's Sapelo high-performance computing cluster due to the massive CPU time required for the permutations)
# outprefix=$gwasdir/4u_gwas
# gwas_blups=$gwasdir/4u_combined_good_traits.update_flowering.tassel.txt
# flowering_covariate=$gwasdir/4u_flowering_covariate.interpolated.tassel.txt
# $scriptdir/4u_RunTasselGwasByChrom.sh $scriptdir $gwas_blups $outdir/chr $outdir $outprefix $flowering_covariate flowering_time $chromlengths $ncores "$TASSEL5" $permdir $empirical_pval_cutoff




# After the above, we have GWAS results on all the original data plus 100 permuted genotype GWAS for each chromosome, yielding empirical p-values
# Since many SNPs are redundant, cluster them according to the methodology of geneSLOPE (Brzyski et al 2015, Controlling the rate of GWAS false discoveries, Genetics 205(1):61-75, doi:  10.1534/genetics.116.193987 )
#    (Had to write own clustering to give enough freedom to choose clustering parameters, etc)

# Gather all SNP hits into a single file to reduce computation
outdir=$gwasdir/4u_gwas # Repeated from above (line ~272) just so can comment out all the above code for debugging
Rscript $scriptdir/4w_CombineParsedGwasResults.r -i $outdir/4v_gwas_results/*.p_adjusted.txt -o $gwasdir/4w_combined_gwas
$TASSEL5 -h $outdir/all_genos.hmp.txt.gz -includeSiteNamesInFile $gwasdir/4w_combined_gwas.markers.txt -export $gwasdir/4w_combined_gwas.hits.hmp.txt.gz
$TASSEL5 -h $gwasdir/4w_combined_gwas.hits.hmp.txt.gz -ld -ldType All -ldHetTreatment Homozygous -export $gwasdir/4w_combined_gwas.hits.ld.txt

# Cluster hits within each result to get the most significant SNP within an LD block
Rscript $scriptdir/4x_ClusterSnpHits.r -i $gwasdir/4w_combined_gwas.combined.txt --genos $gwasdir/4w_combined_gwas.hits.hmp.txt.gz --ld $gwasdir/4w_combined_gwas.hits.ld.txt \
  --cutoff $cluster_cutoff -o $gwasdir/4x_gwas_hits.clustered
