#! /usr/bin/Rscript

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="INput TASSELphenotype file")
parser$add_argument("-n", "--num-perms", help="Number of permutations", type="integer", default=100)
parser$add_argument("-s", "--seed", help="Number of permutations", type="integer", default=1)
parser$add_argument("-o", "--outfile")
args=parser$parse_args()
#setwd("/home/jgw87/Desktop/Working_Files/Metagenomics/16sLeafAmplification/20150601_Batch1_Miseq/")
#args=parser$parse_args(c("-i","5c_sample_phenos.txt","-o","99_perms"))

cat("Permuting phenotypes in",args$infile,"\n")
phenos=read.delim(args$infile, skip=2, header=T)
set.seed(args$seed)
cat(colnames(phenos),"\n")
# Go through and permute each phenotype separately
for(i in 2:ncol(phenos)){
	mypheno=phenos[,c(1,i)]
	unpermuted=2	# Which column is the original
	good_phenos = which(!is.na(mypheno[,unpermuted]))	# Only shuffle the non-NA phenotypes
	
	# Make permutations
	for(n in 1:args$num_perms){
		target_col = unpermuted + n
		mypheno = cbind(mypheno, mypheno[,unpermuted])
		names(mypheno)[target_col] = paste(names(mypheno)[unpermuted], "_perm",n,sep="")
		mypheno[good_phenos,target_col] = mypheno[sample(good_phenos, size=length(good_phenos)),target_col]
	}
	
	# Write out
	write("<Phenotype>",file=args$outfile)
	write(c("taxa", rep("data", ncol(mypheno)-1)), file=args$outfile, ncol=ncol(mypheno), append=T, sep='\t')
	suppressWarnings(write.table(mypheno, file=args$outfile, sep="\t", quote=F, row.names=F, col.names=T, append=T))
# 	break
}
