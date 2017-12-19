#! /bin/bash

# Run Mixed Linear Model GWAS in TASSEL, calculate empirical p-values for significant results, and plot the results for all traits

scriptdir=$1
blups=$2
geno_prefix=$3
targetdir=$4
outprefix=$5
covariate_file=$6
exception=$7 # name of phenotype that is an exception and shouldn't have covariates (flowering time, in my case)
chromlengths=$8
maxprocs=$9
TASSEL5="${10}"
permdir=${11}   # Directory where permuted genotypes are stored
empirical_pval_cutoff=${12}

# Set up directories
blupdir=$targetdir/4u_blups
gwasdir=$targetdir/4u_orig_gwas
permoutdir=$targetdir/4u_permuted_gwas
plotdir=$targetdir/4v_gwas_results
if [ ! -e $blupdir ]; then mkdir $blupdir; fi
if [ ! -e $gwasdir ]; then mkdir $gwasdir; fi
if [ ! -e $permoutdir ]; then mkdir $permoutdir; fi
if [ ! -e $plotdir ]; then mkdir $plotdir; fi


# # Split BLUPs by trait
# nphenos=`head -n3 $blups | tail -n1 | tr '\t' '\n' | wc -l` # Count phenotypes  
# for i in `seq 2 $nphenos`; do
# 
#   # Split BLUPs into individual files
#   name=`head -n3 $blups | tail -n1 | cut -f$i`
#   myblups=$blupdir/4u_blups.$name.txt
#   cut -f1,$i $blups > $blupdir/4u_blups.$name.txt
# 
# #   break
# done


# # # Run GWAS via MLM, including random permutations
# gwas_commands=$outprefix.gwas_commands.txt
# clean_commands=$outprefix.cleanup_commands.txt
# cut_commands=$outprefix.cut_commands.txt
# echo "echo 'Running GWAS by chromosome'" > $gwas_commands
# echo "echo 'Cleaning up MLM output files'" > $clean_commands
# echo "echo 'Cutting out just p-values'" > $cut_commands
# for blups in $blupdir/4u_blups*.txt; do
#   # Extract components. (Could have put this in the earlier loop, but I wanted to keep the tasks separate for clarity)
#   name=${blups/*4u_blups./}
#   name=${name/.txt/}
# 
# #   if [ $name != "flowering_time" ]; then continue; fi   # For Debugging
#   
#   # File names from full OTU pathway cause issues with TASSEL, so cut down
#   if [[ $name == *qiime* ]]; then   
#     name=trait$i;   # NOTE: this is a little unstable if the BLUPS were not split b/c $i isn't defined if that part is commented out
#     i=$(($i+1))
#   fi  
# 
# 
#   # Go through 1 chromosome at a time
#   for chr in $(seq 1 10); do
#     # Files to reference
#     genos=${geno_prefix}$chr.hmp.txt.gz
#     kinship=${geno_prefix}$chr.other_chroms.kinship.txt
#     prefix=$gwasdir/${name}_chr${chr}
#     sitefile=$prefix.sitefile.txt
#     allelefile=$prefix.allelefile.txt
#     pvals=$prefix.pvals.txt
#     
#     # Make gwas commands to be processed in parallel later
#     if [ $name != "$exception" ] ; then 
#       echo "$TASSEL5 -fork1 -h $genos -fork2 -t $blups -fork3 -k $kinship -fork99 -t $covariate_file -combine4 -input1 -input2 -input99 -intersect -combine5 -input3 -input4 \
#         -mlm -mlmMaxP 1 -mlmCompressionLevel None -export ${prefix}_ -runfork1 -runfork2 -runfork3 -runfork99" >> $gwas_commands
#     else
#       echo "$TASSEL5 -fork1 -h $genos -fork2 -t $blups -fork3 -k $kinship                            -combine4 -input1 -input2          -intersect -combine5 -input3 -input4 \
#         -mlm -mlmMaxP 1 -mlmCompressionLevel None -export ${prefix}_ -runfork1 -runfork2 -runfork3" >> $gwas_commands
#     fi  
#     
#     # Clean up output files
#     echo "mv ${prefix}_1.txt $prefix.residuals.txt && mv ${prefix}_2.txt $sitefile && mv ${prefix}_3.txt $allelefile" >> $clean_commands
#     
#     #   Cut out just p-values and gzip other files to save space
#     echo "cut -f1-4,7 $sitefile > $pvals && gzip -f $sitefile && gzip -f $allelefile" >> $cut_commands
#     
#     
#     # Now for permuted values, which are highly similar workflow. I'm sure I could have made this more elegant by calling another script, but I figured the workflow was convoluted enough
#     for perm_genos in $permdir/chr$chr.perm*.hmp.txt.gz; do
#       perm=${perm_genos/*perm/}
#       perm=${perm/.hmp.txt.gz/}
#       
#       # Output files
#       perm_prefix=$permoutdir/${name}_chr${chr}_perm$perm
#       perm_sitefile=$perm_prefix.sitefile.txt
#       perm_allelefile=$perm_prefix.allelefile.txt
#       perm_pvals=$perm_prefix.pvals.txt
#       
#       # For the perms, only keep p-values lower than 0.01 because we're using these as empirical cutoffs; don't need the whole distribution
#       if [ $name != "$exception" ] ; then 
#         echo "$TASSEL5 -fork1 -h $perm_genos -fork2 -t $blups -fork3 -k $kinship -fork99 -t $covariate_file -combine4 -input1 -input2 -input99 -intersect -combine5 -input3 -input4 \
#           -mlm -mlmMaxP 0.01 -mlmCompressionLevel None -export ${perm_prefix}_ -runfork1 -runfork2 -runfork3 -runfork99" >> $gwas_commands
#       else
#         echo "$TASSEL5 -fork1 -h $perm_genos -fork2 -t $blups -fork3 -k $kinship                            -combine4 -input1 -input2          -intersect -combine5 -input3 -input4 \
#           -mlm -mlmMaxP 0.01 -mlmCompressionLevel None -export ${perm_prefix}_ -runfork1 -runfork2 -runfork3" >> $gwas_commands
#       fi  
#       
#       # Clean up output files
#       echo "mv ${perm_prefix}_1.txt $perm_prefix.residuals.txt && mv ${perm_prefix}_2.txt $perm_sitefile && mv ${perm_prefix}_3.txt $perm_allelefile" >> $clean_commands
#       
#       #   Cut out just p-values and gzip other files to save space
#       echo "cut -f1-4,7 $perm_sitefile > $perm_pvals && gzip -f $perm_sitefile && gzip -f $perm_allelefile" >> $cut_commands
#     
# #       break
#       
#     done
#     
#     
#     
# #     break
#   done
# 
# #   break
# done
# cat $gwas_commands | parallel --progress --max-procs $maxprocs 
# cat $clean_commands | parallel --progress --max-procs $maxprocs 
# cat $cut_commands | parallel --progress --max-procs $maxprocs 


# # Gather permuted p-values together and get the most significant p-value per chromosome per trait
# collate_commands=$outprefix.collate_commands.txt
# echo "echo 'Collating p-values from permutations'" > $collate_commands
# for trait in $gwasdir/*_chr1.pvals.txt; do
#   trait=`basename $trait`
#   trait=${trait/_chr1.pvals.txt/}
#   echo $trait
#   echo "Rscript $scriptdir/4u_CollatePermutedPvals.r -i $permoutdir/${trait}_chr*_perm*.pvals.txt -o $plotdir/$trait.pvals.txt" >> $collate_commands
# done
# cat $collate_commands | parallel --progress

# # # # Parse and plot GWAS results with empirical p-values
# plot_commands=$outprefix.plot_commands.txt
# echo "echo Plotting empirical GWAS results" > $plot_commands
# for trait in $gwasdir/*_chr1.pvals.txt; do
#   trait=`basename $trait`
#   trait=${trait/_chr1.pvals.txt/}
#   echo "python3 $scriptdir/4v_ParseGwasResults.py -i $gwasdir/${trait}_chr*.pvals.txt --perms $plotdir/$trait.pvals.txt --sparsify 0.1 --empirical-p-cutoff $empirical_pval_cutoff -o $plotdir/$trait.parsed --chromlengths $chromlengths" >> $plot_commands
# #   break
# done
# cat $plot_commands | parallel --progress