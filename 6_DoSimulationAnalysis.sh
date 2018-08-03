#! /bin/bash

# Run a simulation analysis with different QTL settings to make sure these methods can actually find causal QTN

scriptdir=$1
simdir=$2
genotypes=$3
phenotypes=$4
covariate_file=$5
traits="$6"
num_qtn="$7"
herits="$8"
nreps=$9
gwas_outdir=${10}
genos_bychrom_prefix=${11}
chromlengths=${12}
permdir=${13}
empirical_pval_cutoff=${14}
maxprocs=${15}
TASSEL5="${16}"

splitdir=$simdir/6a_split_genos
myperms=$simdir/6a_permuted_genos
if [ ! -e $splitdir ]; then mkdir $splitdir; fi
if [ ! -e $myperms ]; then mkdir $myperms; fi

# # Pull phenotypes out from GWAS input file
Rscript $scriptdir/6a_SubsetPhenotypes.r -i $phenotypes -t $traits -o $simdir/6a_target_phenos

# Get BLUPs and residuals after fitting a kinship matrix
$TASSEL5 -h $genotypes -ck -export $simdir/6b_kinship.txt
zcat $genotypes | head -n 2 > $simdir/6b_dummy_genos.hmp.txt

# Subset to just the genotypes that had pvals of 0.001 or better among all real phenotypes. (Meant to speed computation; presumably other sites are mostly irrelevant/underpowered)
python3 $scriptdir/6a_IdentifyMeaningfulSites.py -i $gwas_outdir/*.pvals.txt --max-p 0.001 -o $simdir/6a_target_sites.txt
$TASSEL5 -h $genotypes -includeSiteNamesInFile $simdir/6a_target_sites.txt -export $simdir/6a_genos_short.hmp.txt.gz

# Separate out the above smaller genotype files into individual chromosomes, but still keep the original kinship matrices (since filtering genotypes is just meant to have MLM take less time)
for chr in `seq 1 10`; do
    cp ${genos_bychrom_prefix}${chr}.other_chroms.kinship.txt $splitdir # Copy over original kinship matrix of all other chromosomes
    $TASSEL5 -h $simdir/6a_genos_short.hmp.txt.gz -separate $chr -export $splitdir/chr$chr.hmp.txt.gz   # Separate out per-chromosome genotypes; file name should basically match that of the genos_bychrom_prefix, just in a different directory
done

# Filter all the permuted genotype files to just the above sites
for hmp in $permdir/*.hmp.txt.gz; do
    filename=`basename $hmp`
    $TASSEL5 -h $hmp -includeSiteNamesInFile $simdir/6a_target_sites.txt -export $myperms/$filename
done

# Make numeric genotypes to do simulations from
$TASSEL5 -h  $simdir/6a_genos_short.hmp.txt.gz -NumericalGenotypePlugin -endPlugin -ExportPlugin -saveAs $simdir/6b_genos.numeric.txt -format ReferenceProbability

# Run analysis separately for each traits
seed=0  # Random seed for permuting residuals
for trait in $traits; do
    seed=$(($seed + 1))
    genos=$simdir/6b_dummy_genos.hmp.txt    # Single site since only need residuals
    kinship=$simdir/6b_kinship.txt
    blups=$simdir/6a_target_phenos.$trait.tassel.txt

    # Create subdirectory
    workdir=$simdir/6c_$trait
    if [ ! -e $workdir ]; then mkdir $workdir; fi
    
    # Run MLM to get residuals
    if [ $trait != "flowering_time" ] ; then 
      $TASSEL5 -fork1 -h $genos -fork2 -t $blups -fork3 -k $kinship -fork99 -t $covariate_file -combine4 -input1 -input2 -input99 -intersect -combine5 -input3 -input4 \
        -mlm -mlmMaxP 1 -mlmCompressionLevel None -export $workdir/6c_mlm -runfork1 -runfork2 -runfork3 -runfork99
    else
      $TASSEL5 -fork1 -h $genos -fork2 -t $blups -fork3 -k $kinship                            -combine4 -input1 -input2          -intersect -combine5 -input3 -input4 \
        -mlm -mlmMaxP 1 -mlmCompressionLevel None -export $workdir/6c_mlm -runfork1 -runfork2 -runfork3
    fi 
    
    # Rename residuals
    mv $workdir/6c_mlm1.txt $workdir/6c_mlm.residuals.txt
    mv $workdir/6c_mlm2.txt $workdir/6c_mlm.stats.txt
    
    # Scramble residuals
    Rscript $scriptdir/6d_PermuteResiduals.r -i $workdir/6c_mlm.residuals.txt -o $workdir/6d_mlm.residuals.permuted.txt --seed $seed --reps $nreps
    
    # Add QTL
    Rscript $scriptdir/6e_SimulateQtl.r -i $workdir/6d_mlm.residuals.permuted.txt --genos $simdir/6b_genos.numeric.txt -o $workdir/6e_residuals_with_qtl --num-qtl $num_qtn --herits $herits --seed $seed
    
    # Add these new residuals back in place of the old ones to get phenotypes ready for complete GWAS analysis
    Rscript $scriptdir/6f_RegenerateBlups.r --raw $blups --residuals $workdir/6c_mlm.residuals.txt --perms $workdir/6e_residuals_with_qtl.phenos.txt -o $workdir/6f_blups_with_qtl.txt
    
#     break
done


# Now go through and run GWAS (Put in a separate loop for ease of debugging and analysis
for trait in $traits; do
    workdir=$simdir/6c_$trait
    genos=$simdir/6b_dummy_genos.hmp.txt
    kinship=$simdir/6b_kinship.txt
    blups=$workdir/6f_blups_with_qtl.txt

    extradir=$workdir/6g_extra_files
    if [ ! -e $extradir ]; then mkdir $extradir; fi
    
    # Run MLM on simulated data to get heritabilities
    if [ $trait != "flowering_time" ] ; then 
      $TASSEL5 -fork1 -h $genos -fork2 -t $blups -fork3 -k $kinship -fork99 -t $covariate_file -combine4 -input1 -input2 -input99 -intersect -combine5 -input3 -input4 \
        -mlm -mlmMaxP 1 -mlmCompressionLevel None -export $workdir/6g_mlm -runfork1 -runfork2 -runfork3 -runfork99
    else
      $TASSEL5 -fork1 -h $genos -fork2 -t $blups -fork3 -k $kinship                            -combine4 -input1 -input2          -intersect -combine5 -input3 -input4 \
        -mlm -mlmMaxP 1 -mlmCompressionLevel None -export $workdir/6g_mlm -runfork1 -runfork2 -runfork3
    fi 
    
    # Parse output files (I really wish TASSEL let me set them individually)
    for outfile in $workdir/6g_mlm*.txt; do
        is_stats=`head -n 1 $outfile | grep "MarkerR2" | wc -l | cut -f1`
        if [[ $is_stats -eq 1 ]] ; then 
            mv $outfile $outfile.stats
        else
            mv $outfile $extradir
        fi
    done
    
    # Extract heritabilities to compare original with new
    python3 $scriptdir/6h_PlotPermHeritabilities.py --orig $workdir/6c_mlm.stats.txt --perms $workdir/6g_mlm*.txt.stats -o $workdir/6h_perm_heritabilities.png
    
#     # Run on full dataset and plot GWAS
    outprefix=$workdir/6i_gwas
    genos_prefix=$splitdir/chr
    outdir=$outprefix
    if [ ! -e $outdir ]; then mkdir $outdir; fi
    $scriptdir/4u_RunTasselGwasByChrom.sh $scriptdir $blups $genos_prefix $outdir $outprefix $covariate_file flowering_time $chromlengths $maxprocs "$TASSEL5" $myperms $empirical_pval_cutoff
    
    # Subset original hapmap to just simulated QTL so have their physical positions
    Rscript -e "x=read.delim('$workdir/6e_residuals_with_qtl.simulation_stats.txt', stringsAsFactors=F); qtn=x[,'qtn']; qtn=unlist(lapply(qtn, strsplit, split=',')); write(sort(unique(qtn)), '$workdir/6j_qtn.txt')"
    $TASSEL5 -h $genotypes -includeSiteNamesInFile $workdir/6j_qtn.txt -export $workdir/6j_qtn.hmp.txt
    
    # Replot Manhatten plots with simulated QTL marked (code modified from 4u_RunTasselGwasByChrom.sh; overwrites ones made by the above script (but are identical except for QTL locations anyway)
    plot_commands=$outprefix.plot_with_qtl_commands.txt
    gwasdir=$outdir/4u_orig_gwas
    plotdir=$outdir/4v_gwas_results
    echo "echo Plotting empirical GWAS results with QTL marked" > $plot_commands
    for trait in $gwasdir/*_chr1.pvals.txt; do
        trait=`basename $trait`
        trait=${trait/_chr1.pvals.txt/}
        echo "python3 $scriptdir/6j_ParseGwasResultsWithSimulatedQtl.py -i $gwasdir/${trait}_chr*.pvals.txt --perms $plotdir/$trait.pvals.txt --sparsify 0.1 --empirical-p-cutoff $empirical_pval_cutoff -o $plotdir/$trait.parsed_with_qtl --chromlengths $chromlengths --qtl $workdir/6e_residuals_with_qtl.simulation_stats.txt --hapmap $workdir/6j_qtn.hmp.txt" >> $plot_commands
    #   break
    done
    cat $plot_commands | parallel --progress

    # Parse QTL results down to a single file
    python3 $scriptdir/6k_CollateQtlResults.py -i $plotdir/*.parsed_with_qtl.found_qtl.txt --qtl $workdir/6e_residuals_with_qtl.simulation_stats.txt -o $workdir/6k_qtl_simulation_results -p 0.1 0.05 0.01
    
    
#     break  
done
