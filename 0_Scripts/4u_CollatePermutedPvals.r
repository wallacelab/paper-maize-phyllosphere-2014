#! /usr/bin/Rscript
 
library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infiles", nargs="*", help="List of files containing significant p-values from permutations; should be format '*.permX.pvals.txt', where X is an integer")
parser$add_argument("-o", "--outfile")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/282MaizeLeaf16s_cleaned/4_GWAS/4u_gwas/4u_permuted_gwas/COG0303_2759/")
# args=parser$parse_args(c("-i","COG0303_2759.chr1.perm0.pvals.txt", "COG0303_2759.chr1.perm1.pvals.txt", "COG0303_2759.chr1.perm2.pvals.txt", "COG0303_2759.chr2.perm2.pvals.txt", "-o","99_tmp.txt"))

# Load data
cat("Collating p-values from",length(args$infiles),"input files\n")
pvals = lapply(args$infiles, read.delim)

# Cut down to most significant
pvals = lapply(pvals, function(x){
  x$Pos=NULL	# Drop "Pos" column
  x=x[which.min(x$p),]
  return(x)
})
pvals=do.call(rbind, pvals)

# Make output by putting each chromosome in a colunm, padding with NAs if necessary
output = split(pvals, pvals$Chr)
max_length = max(unlist(sapply(output, nrow)))
output = lapply(output, function(x){
  return(x$p[1:max_length])
})
output=do.call(cbind, output)
colnames(output) = paste("chr",colnames(output), sep="")

# Write out
cat("\tWriting output of",nrow(output),"permutations across",ncol(output),"chromosomes to",args$outfile,"\n")
write.table(output, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T)