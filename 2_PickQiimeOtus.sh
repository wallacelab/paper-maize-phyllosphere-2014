#! /bin/bash

# Run QIIME on the maize microbiome data

scriptdir=$1
qiimedir=$2
parsedir=$3
plate_key=$4
read_dist=$5
min_otu_size=$6
sample_threshold=$7
ignore_samples=$8
maxprocs=$9
KRONA=${10}

###############
# Run QIIME pipeline
###############

pickdir=$qiimedir/2c_PickOtus

# Make QIIME mapping key; note that the inbred pop key is mostly based off Flint-Garcia 2005's supplemental table, but with some names changed to match and a few added that weren't in the original
Rscript $scriptdir/2a_MakeQiimeKey.r -i $plate_key -m $qiimedir/2a_qiime_sample_key.tsv -l $scriptdir/0_inbred_subpops.txt #--raw-depth-only -r $read_dist
validate_mapping_file.py -m $qiimedir/2a_qiime_sample_key.tsv --disable_primer_check --not_barcoded --added_demultiplex_field=SampleID
mv -f 2a_qiime_sample_key.tsv* $qiimedir   # Move the validation files that are created in the root directory to the qiime directory

# # Change underscores to periods in sample names because QIIME doesn't accept them for some reason. Also add read numbers
python3 $scriptdir/2b_FixSampleNames.py -i $parsedir/*/joined_with_sample_name.fasta.gz  -o $qiimedir/2b_seqs_combined.fixed.fna #--debug

# OTU picking; parameter file just has a flag to enable reverse strand matching
# Would be nice to do this in parallel, but QIIME keeps failing to execute some threads (I don't know why)
pick_open_reference_otus.py -o $pickdir -i $qiimedir/2b_seqs_combined.fixed.fna -p $scriptdir/2c_otu_picking_params.txt --suppress_step4 -f --min_otu_size $min_otu_size

# Pynast alignment keeps messing up one clade of methylobacteria, so use ClustalOmega to get a more accurate alignment
clustalo -i $pickdir/rep_set.fna -o $pickdir/rep_set_aligned.clustalo.fasta --threads=$maxprocs --verbose
make_phylogeny.py -i $pickdir/rep_set_aligned.clustalo.fasta -o $pickdir/rep_set.clustalo.tre

# Identify good samples to keep / bad samples to remove
ignore_samples=`echo $ignore_samples| tr '_' '.'` # Again, change underscores to periods
biom summarize-table -i $pickdir/otu_table_mc${min_otu_size}_w_tax_no_pynast_failures.biom > $qiimedir/2c_biom_summary.samples.txt
python3 $scriptdir/2d_IdentifyGoodTaxaInBiomFile.py -i $qiimedir/2c_biom_summary.samples.txt --blank-patterns "blank" "BLANK" --sample-to-blank-threshhold $sample_threshold \
  --remove-samples $ignore_samples -o $qiimedir/2d_otu_table_samples_to_keep.txt -f $qiimedir/2d_otu_table_samples_to_remove.txt
filter_samples_from_otu_table.py -i $pickdir/otu_table_mc${min_otu_size}_w_tax_no_pynast_failures.biom -o $qiimedir/2e_otu_table.sample_filtered.biom --sample_id_fp $qiimedir/2d_otu_table_samples_to_keep.txt

# Remove mitochonrdia and chloroplast (shouldn't be many of latter due to using anti-chloroplast primers)
filter_taxa_from_otu_table.py -i $qiimedir/2e_otu_table.sample_filtered.biom -o $qiimedir/2f_otu_table.sample_filtered.no_mitochondria_chloroplast.biom --negative_taxa f__mitochondria,c__Chloroplast

# Make summaries and a text version
prefix=$qiimedir/2f_otu_table.sample_filtered.no_mitochondria_chloroplast
biom summarize-table -i $prefix.biom -o $prefix.samples.txt
biom summarize-table -i $prefix.biom --observations -o $prefix.observations.txt
biom convert -i $prefix.biom -o $prefix.txt --to-tsv

# Copy accessory files to same directory and with same prefix for easier access
cp $pickdir/rep_set.clustalo.tre $prefix.tre
cp $pickdir/uclust_assigned_taxonomy/rep_set_tax_assignments.txt $prefix.taxonomy.txt
cp $qiimedir/2a_qiime_sample_key.tsv $prefix.key.tsv

# Gzip the biggest files to save space
gzip $qiimedir/2b_seqs_combined.fixed.fna $pickdir/final_otu_map* $pickdir/*/failures.fasta $pickdir/*/failures*.uc

# Make Krona plots of the output
python3 $scriptdir/2g_FormatTaxonomyForKrona.py -i $prefix.taxonomy.txt -b $prefix.biom  -o $qiimedir/2g_krona_input.txt --include-otus
$KRONA -o $qiimedir/2g_krona_plot.html $qiimedir/2g_krona_input.txt