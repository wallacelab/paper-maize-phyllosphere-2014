#! /bin/bash

# Gather maize genotype data for use in GWAS analysis
# NOTE: This is mostly included for reference; the filtered genotype data should be included as part of the pipeline

scriptdir=$1
genodir=$2
plate_key=$3
source_genos=$4
source_keyfile=$5
additional_samples=$6
synonyms=$7
TASSEL5=$8
min_freq=$9
min_count=${10}


# Make keyfile of target lines and subset out just those taxa; a bit complicated b/c the full hapmap file won't fit in memory
taxalist=$genodir/0e_samples_gbs_taxa.txt
grep -e '282 maize diversity panel' -e 'Ames282' $source_keyfile | cut -f1 > $taxalist  # These are how this population is labeled in the Panzea GBS keyfile
cut -f1 $source_keyfile | grep -i --file=$additional_samples  >> $taxalist   # For the ones that aren't part of the 2 groups listed above
python3 $scriptdir/0e_SubsetHapmapByTaxa.py -i $source_genos -o $genodir/0e_samples.hmp.txt.gz -t $taxalist # --debug
$TASSEL5 -SortGenotypeFilePlugin -inputFile $genodir/0e_samples.hmp.txt.gz  -outputFile  $genodir/0f_samples.sorted.hmp.txt.gz #TASSEL is picky about order, and the original file is in a different order than it wants

# Filter genotypes to remove duplicate taxa (taking the ones with highest coverage, or--failing that--the ones with less missing data) 
$TASSEL5 -h $genodir/0f_samples.sorted.hmp.txt.gz -genotypeSummary taxa -export $genodir/0f_samples.sorted.taxasummary.txt  
python3 $scriptdir/0g_SelectTaxaOnCoverage.py -t $genodir/0f_samples.sorted.taxasummary.txt -o $genodir/0g_samples_best_taxa.txt --best $genodir/0g_samples_best_taxa.names.txt

# Filter genotypes to remove low-coverage sites
$TASSEL5 -h $genodir/0f_samples.sorted.hmp.txt.gz -includeTaxaInFile $genodir/0g_samples_best_taxa.names.txt -filterAlign -filterAlignMinFreq $min_freq -filterAlignMinCount $min_count \
  -export $genodir/0h_samples_sorted_filtered.hmp.txt.gz

# Make kinship matrix
$TASSEL5 -h $genodir/0h_samples_sorted_filtered.hmp.txt.gz -ck -export $genodir/0h_samples_sorted_filtered.kinship.txt

