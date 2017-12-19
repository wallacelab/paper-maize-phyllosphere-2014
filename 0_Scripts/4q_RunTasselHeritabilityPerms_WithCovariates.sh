#! /bin/bash

# Run TASSEL heritability permutations without any covariates

scriptdir=$1
blups=$2
dummy_genos=$3
kinship=$4
nperms=$5
cutoff=$6
targetdir=$7
outprefix=$8
covariate_file=$9
exception=${10} # name of phenotype that is an exception and shouldn't have covariates
TASSEL5="${11}"


parallel_mlm=""
parallel_perms=""
nphenos=`head -n3 $blups | tail -n1 | tr '\t' '\n' | wc -l` # Count phenotypes
for i in `seq 2 $nphenos`; do
  # Setup
  name=`head -n3 $blups | tail -n1 | cut -f$i`
  cut -f1,$i $blups > $targetdir/4q_blups.$name.txt
  outdir=$targetdir/$name
  if [ ! -e $outdir ] ; then mkdir $outdir; fi
  
  # Load permutation commands for parallel processing later
  parallel_perms="$parallel_perms \n Rscript $scriptdir/4q_PermutePhenotypesForHeritability.r -i $targetdir/4q_blups.$name.txt -o $targetdir/4q_blups.$name.permuted.txt -n $nperms"
   
  # Load MLM commands for parallel processing later
  statsfile=$(($nperms+2))
  if [ $name != "$exception" ] ; then 
    parallel_mlm="$parallel_mlm \n $TASSEL5 -fork1 -h $dummy_genos -fork2 -t $targetdir/4q_blups.$name.permuted.txt -fork3 -k $kinship -fork99 -t $covariate_file -combine4 -input1 -input2 -input99 -intersect -combine5 -input3 -input4 \
      -mlm -mlmMaxP 1 -mlmCompressionLevel None -export $outdir/4q_$name -runfork1 -runfork2 -runfork3 -runfork99; mv $outdir/*${statsfile}.txt $outdir/4r_${name}_stats.txt"    # No compression because causes TASSEL to hang on a lot of the permutations
  else
    parallel_mlm="$parallel_mlm \n $TASSEL5 -fork1 -h $dummy_genos -fork2 -t $targetdir/4q_blups.$name.permuted.txt -fork3 -k $kinship -combine4 -input1 -input2 -intersect -combine5 -input3 -input4 \
      -mlm -mlmMaxP 1 -mlmCompressionLevel None -export $outdir/4q_$name -runfork1 -runfork2 -runfork3; mv $outdir/*${statsfile}.txt $outdir/4r_${name}_stats.txt"
  fi  
  
done
# Run stored parallel commands
echo "Running stored MLM commands in parallel for $blups"
echo -e $parallel_perms | parallel --progress
echo -e $parallel_mlm   | parallel --progress
# 
# Remove residual files to not clutter up computer
echo "Removing residuals files from subdirectories"
for subdir in $targetdir/*/; do
  rm $subdir/4q_*.txt
done
  
# Plot heritability distributions and filter for ones with good empirical p-values
echo "Plotting heritability distributions"
python3 $scriptdir/4r_CalculateHeritability.py -i $targetdir/*/4r_*_stats.txt -o $outprefix.txt -g $outprefix.png --passfile $outprefix.good_traits.txt --pval-cutoff $cutoff


