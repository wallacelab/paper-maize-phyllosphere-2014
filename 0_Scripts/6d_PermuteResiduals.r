#! /usr/bin/env Rscript

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infile")
parser$add_argument("-o", "--outfile")
parser$add_argument("-r", "--reps", type="integer", default=3, help="Number of independent permutations to run")
parser$add_argument("-s", "--seed", type="integer", default=1, help="Random seed for permutations")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/282MaizeLeaf16s_cleaned/6_PowerSimulation/6c_flowering_time/")
# args=parser$parse_args(c("-i","6c_mlm.residuals.txt","-o","99_tmp"))

# Load data
cat("Permuting residuals in", args$infile,"\n")
data = read.delim(args$infile, skip=2, header=T, check.names=F)
if(ncol(data) != 2){stop("This script assumes a single phenotype per TASSEL-formatted input file\n")}

# Permute data
set.seed(args$seed)
perms = sapply(1:args$reps, function(i){
    sample(data[,2], size=nrow(data), replace=F)
})
colnames(perms) = paste(names(data)[2], "_perm", 1:ncol(perms), sep="")

# Create output
output = data.frame(Taxa=data$Taxa, perms)

# Write
datatypes = paste("taxa", paste(rep("data", ncol(perms)), collapse="\t"), sep="\t")
write("<Phenotype>", args$outfile)
write(datatypes, args$outfile, append=T)
write.table(output, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T, append=T)
