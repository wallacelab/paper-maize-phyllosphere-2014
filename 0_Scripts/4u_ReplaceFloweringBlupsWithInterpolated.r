#! /usr/bin/Rscript
 
library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--interpfile", help="TASSEL covariate file with interpolated flowering time BLUPs")
parser$add_argument("-b", "--blupfile", help="BLUP file to update")
parser$add_argument("-o", "--outfile", help="Updated BLUP file")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/282MaizeLeaf16s_cleaned/4_GWAS/")
# args=parser$parse_args(c("-i","4u_flowering_covariate.interpolated.tassel.txt", "-b", "4t_combined_good_traits.tassel.txt", "-o","99_tmp.txt"))

# Load data
cat("Updating BLUP values in",args$blupfile,"with interpolated flowering time from",args$interpfile,"\n")
blup_header=scan(args$blupfile, what=character(), sep='\n', nlines=2)
blups=read.delim(args$blupfile, skip=2, header=T, check.names=F)
flowering=read.delim(args$interpfile, skip=2, header=T, check.names=F)

# Match up (should be already, but to be sure)
mymatch = match(blups$Taxa, flowering$Taxon)
if(!identical(blups$Taxa[mymatch], flowering$Taxon)){
  cat("\tWARNING! Mismatches found when matching taxa names!\n")
}

# Swap over, saving old values for some sanity checks
old = blups$flowering_time
blups$flowering_time = flowering$flowering_time[mymatch]

# Sanity checks
corr = cor(old, blups$flowering_time, use='pairwise')
before=sum(is.na(old))
after=sum(is.na(blups$flowering_time))
cat("\tCorrelation between new and old values is",corr,"\n")
cat("\tBefore interpolation had",before,"NAs; afterward have",after,"\n")

# Write new blups out
cat("Writing to new BLUP file (this will complain about adding column headers)\n")
write(blup_header, file=args$outfile)
write.table(blups, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T, append=T)
