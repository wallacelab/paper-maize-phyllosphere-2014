#! /usr/bin/Rscript

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infiles", nargs="*", help="Input files (results of 4v_ParseGwasResults.py, with the '.p_adjusted.txt' suffix")
parser$add_argument("-o", "--outprefix")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/282MaizeLeaf16s_cleaned/4_GWAS/")
# args=parser$parse_args(c("-i","4u_gwas/4v_gwas_results/COG0067_1122.parsed.p_adjusted.txt","4u_gwas/4v_gwas_results/COG0070_785.parsed.p_adjusted.txt", "4u_gwas/4v_gwas_results/COG0784_2804.parsed.p_adjusted.txt", "-o","99_tmp"))

# Load data
cat("Loading parsed GWAS results from",length(args$infiles),"input files\n")
input = lapply(args$infiles, read.delim, row.names=1, stringsAsFactors=F)

# Output combined data frame
combined_out = paste(args$outprefix, ".combined.txt", sep="")
cat("\tWriting combined output table to",combined_out,"\n")
combined = do.call(rbind, input)
combined=combined[order(combined$Trait, combined$Chr, combined$Pos),]	# sort
write.table(combined, file=combined_out, sep='\t', quote=F, row.names=F, col.names=T)

# Output just SNP names
marker_out = paste(args$outprefix, ".markers.txt", sep="")
cat("\tWriting marker names to",marker_out,"\n")
markers = sort(unique(combined$Marker))
write(markers, file=marker_out)


# Output table of counts for how often specific SNPs were found
count_out = paste(args$outprefix, ".marker_counts.txt", sep="")
cat("\tWriting count of hits per marker to",count_out,"\n")
counts = as.data.frame(table(combined$Marker))
names(counts)=c('Marker', 'Count')
write.table(counts, file=count_out, sep='\t', quote=F, row.names=F, col.names=T)


# Output table of cutoffs for each trait (most involved)
# # First, filter for traits that had at least 1 hit
goodhits = split(combined, combined$Trait)	# Split 'combined' instead of going back to 'input' because it automatically removes any traits with 0 hits
traits=sapply(goodhits, function(x){unique(x$Trait)})
summary=data.frame(traits=traits)

# # Now get the count for each trait at different cutoffs
under_cutoff=function(x, cutoff=1){	# Helper function to determine how many empirical p-values are under a certain cutoff
	return(sum(x$empirical_pval <= cutoff))
}
summary$"p≤0.10" = sapply(goodhits, under_cutoff, cutoff=0.10)
summary$"p≤0.05" = sapply(goodhits, under_cutoff, cutoff=0.05)
summary$"p≤0.01" = sapply(goodhits, under_cutoff, cutoff=0.01)

# # Finally, write out
summary_out = paste(args$outprefix, ".trait_summary.txt", sep="")
cat("\tWriting summary of hits per trait to",summary_out,"\n")
write.table(summary, file=summary_out, sep='\t', quote=F, row.names=F, col.names=T)