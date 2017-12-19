#! /usr/bin/env bash

# Filter a MED abundnace table by both samples and MED nodes and then get diversity statistics on the result (including the QIIME pipeline)
scriptdir=$1
divdir=$2
filestem=$3 # Base file stem for biom file, key, taxonomy, etc.
rarefaction_level=$4   # What to rarefy at for community analyses
core_fraction=$5
KRONA=$6
maxprocs=$7

# Standard QIIME alpha and beta diversity analysis; no parallelization b/c keeps failing to launch some processes
core_diversity_analyses.py -i $filestem.biom -o $divdir/3a_diversity_unifrac -m $filestem.key.tsv --sampling_depth $rarefaction_level --tree_fp $filestem.tre --recover_from_failure 
core_diversity_analyses.py -i $filestem.biom -o $divdir/3a_diversity_bray_curtis -m $filestem.key.tsv --sampling_depth $rarefaction_level --nonphylogenetic  --recover_from_failure

# Identify core microbiome (Note: this is a first-pass approximation. Use the one in 9_PrettifyGraphics for the paper)
nsamples=`grep -c --no-filename "^LMA" $filestem.samples.txt`   # Count number of samples
min_samples=`echo -e $core_fraction \* $nsamples / 1| bc`   # Use bc to do the math because bash can't
filter_otus_from_otu_table.py -i $filestem.biom -o $divdir/3b_core_community.biom --min_samples $min_samples
biom summarize-table -i $divdir/3b_core_community.biom -o $divdir/3b_core_community.otus.txt --observations
summarize_taxa_through_plots.py -i $divdir/3b_core_community.biom -o $divdir/3b_core_taxa_summaries --force

# Make Krona plot of core microbiome
python3 $scriptdir/2g_FormatTaxonomyForKrona.py -i $filestem.taxonomy.txt -b $divdir/3b_core_community.biom -o $divdir/3c_core_community.krona_input.txt --include-otus
$KRONA -o $divdir/3c_core_community.krona.html $divdir/3c_core_community.krona_input.txt
