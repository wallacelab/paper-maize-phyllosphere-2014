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


