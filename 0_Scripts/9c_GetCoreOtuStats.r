#! /usr/bin/env Rscript

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input OTU table (text)")
parser$add_argument("-o", "--outfile", help="Output table of core OTUs")
parser$add_argument("-a", "--allfile", help="Output table of all OTU stats")
parser$add_argument("-p", "--min-presence-across-samples", type="double", default=0.8, help="What fraction of samples an OTU has to be in to be considered 'core'")
parser$add_argument("-c", "--min-count-in-samples", type="integer", default=1, help="How many times an OTU has to be in a given sample to count as present")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/282MaizeLeaf16s_cleaned')
# args=parser$parse_args(c("-i","4_GWAS//4f_orig_otus.filtered.txt","-o",'99_tmp_core.txt','-a', '99_tmp_all.txt'))

# Load data
biom = as.matrix(read.delim(args$infile, skip=1, row.names=1))
cat("Summarizing",nrow(biom),"OTUs across",ncol(biom),"samples\n")

# Make summary data frame
summary=data.frame(otu=rownames(biom))

# Fraction presence across samples
presence = biom >= args$min_count_in_samples
summary$presence = rowSums(presence) / ncol(presence)

# Mean fraction of total sample
totals = matrix(colSums(biom), nrow=nrow(biom), ncol=ncol(biom), byrow=T)
fractions = biom / totals	# Fraction within each sample
summary$mean_fraction_of_sample = rowMeans(fractions)
summary$median_fraction_of_sample = apply(fractions, MARGIN=1, FUN=median)

# Fraction of total reads
summary$fraction_total = rowSums(biom) / sum(biom)

#Write entire OTU summary out
if(!is.null(args$allfile)){
  cat("\tWriting all OTU summaries to",args$allfile,"\n")
  write.table(summary, file=args$allfile, sep='\t', quote=F, row.names=F, col.names=T)
}

# Filter to core OTUs
filtered = subset(summary, summary$presence >= args$min_presence_across_samples)
cat("\tIdentified",nrow(filtered),"core OTUs with presence in at least",args$min_presence_across_samples,"of samples\n")
cat("\t\tThese OTUs are",sum(filtered$fraction_total),"of the total input dataset\n")

# Write out core OTU summary
if(!is.null(args$outfile)){
  cat("\tWriting core OTU summaries to",args$outfile,"\n")
  filtered = filtered[order(filtered$fraction_total, decreasing=T),]
  write.table(filtered, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T)
}