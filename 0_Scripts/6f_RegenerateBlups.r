#! /usr/bin/env Rscript

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--raw", help="Original BLUPs file (in TASSEL format)")
parser$add_argument("-p", "--perms", help="Permuted residuals with QTL to add back")
parser$add_argument("-r", "--residuals", help="Original residuals from mixed linear model")
parser$add_argument("-o", "--outfile", help="Output file")
args=parser$parse_args()
# setwd("/home/jgwall/Projects/282MaizeLeaf16s_cleaned/6_PowerSimulation/6c_flowering_time/")
# args=parser$parse_args(c("-i","../6a_target_phenos.flowering_time.tassel.txt","-r","6c_mlm.residuals.txt", "-p", "6e_residuals_with_qtl.txt", "-o", "99_tmp.txt" ))

# Load data
cat("Loading data to created BLUPs with simulated QTL\n")
blups=read.delim(args$raw, skip=2, row.names=1, header=T)
residuals=read.delim(args$residuals, skip=2, row.names=1, header=T)
perms=read.delim(args$perms, skip=2, row.names=1, header=T)
cat("\tLoaded", ncol(perms),"permuted datasets to work with\n")

# Sanity checks
if(ncol(blups)!=1){stop("This script assumes a single phenotpye in the BLUP file but ",ncol(blups)," found in",args$raw,"\n")}
if(ncol(residuals)!=1){stop("This script assumes a single data column in the original residuals file but ",ncol(residuals)," found in",args$residuals,"\n")}

# Make sure taxa are all in right order
cat("Making sure row names are all matched\n")
taxa=intersect(rownames(blups), intersect(rownames(residuals), rownames(perms)))
blups=blups[match(taxa, rownames(blups)),, drop=F]
residuals = residuals[match(taxa, rownames(residuals)),, drop=F]
perms=perms[match(taxa, rownames(perms)),, drop=F]
# Confirm match worked
if(!identical(rownames(blups), rownames(residuals)) || !identical(rownames(blups), rownames(perms))){
    stop("\t\tUnable to match row names; please check code\n")
}

# Get non-residual parts and add to permuted residuals
modelvals = blups - residuals
newblups = unlist(modelvals) + perms    # Have to unlist so it's a vector; otherwise R complains

# Write out
cat("Writing output to",args$outfile,"\n")
datatypes = paste("taxa", paste(rep("data", ncol(newblups)), collapse="\t"), sep="\t")
output=data.frame(Taxa=rownames(newblups), newblups) # Add header for Taxa column
names(output) = gsub(names(output), pattern=".", repl="_", fixed=T) # Replace periods with underscores so TASSEL handles them properly
write("<Phenotype>", args$outfile)
write(datatypes, args$outfile, append=T)
write.table(output, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T, append=T)

