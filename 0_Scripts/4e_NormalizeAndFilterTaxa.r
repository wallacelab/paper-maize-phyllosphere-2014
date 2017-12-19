#! /usr/bin/Rscript

# Convert OTU counts (transformed or otherwise) into BLUPs

library(argparse)
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input file of traits as a matrix")
parser$add_argument("-o", "--outfile", help="Output file of traits as a matrix")
parser$add_argument("-p", "--printfile", help="File to print diagnostics")
parser$add_argument("--sample-summary", help="Summary of biom file across samples, from BIOM command line. Used to get normalization for each sample")
parser$add_argument("-s", "--skip", type="integer", default=0, help="Lines to skip in the BIOM file (useful if has additional header lines)")
parser$add_argument("-c", "--normalize-count", type="integer", default=1000, help="Count to normalize to")
parser$add_argument("-z", "--min-count-across-samples", type="double", default=0.0, help="Maximum fraction of missing counts (=0) acoss samples for an OTU to be kept, ie, 0.8 means anything with 20% or more 0 counts is filtered out")
parser$add_argument("--skip-normalization", default=FALSE, action="store_true", help="Whether to skip the normalization step")
args=parser$parse_args()
# setwd('/home/jgwall/Desktop/jason/Projects/282MaizeLeaf16s_redux/4_GWAS/')
# args=parser$parse_args(c("-i","4e_combined.collapsed_taxa.txt","-o",'99_tmp.txt','-k','4f_merged_keyfile.qiime.tsv','-n', 'raw_depth', "-p", "99_tmp.png"))
# setwd('/home/jgwall/Projects/282MaizeLeaf16s_redux/4_GWAS_OldDataSummaryStats/')
# args=parser$parse_args(c("-i","4b_qiime.no_mitochondria.filtered20.txt","-o",'99_tmp.txt','-k','4b_qiime.no_mitochondria.filtered20.keyfile.txt','-n', 'nomito_count', "-p", "99_tmp.png","--skip", "1"))

cat('Normalizing taxa in',args$infile,"to",args$normalize_count,"reads\n")
data=read.delim(args$infile, check.names=F, row.names=1, skip=args$skip)
summary=read.table(args$sample_summary, skip = grep("Counts/sample detail", readLines(args$sample_summary)), col.names=c('sample','count'))	# Skip down to the part with count data
summary$sample=sub(summary$sample, pattern=":", repl="")

# Build data frame with the appropriate setup
keymatch = match(colnames(data), summary$sample)
newkey = summary[keymatch,]

# Make matrix of total counts and use to normalize
if(args$skip_normalization){
  cat("NOTE: Sample normalization skipped due to the --skip-normalization flag!\n")
  data.norm=data
} else{
  totals = matrix(newkey$count, nrow=nrow(data), ncol=ncol(data), byrow=T)
  data.norm = data / totals * args$normalize_count
  data.norm=round(data.norm, digits=0)
}

# Filter taxa based on minimum count (same as normalize count)
cat("Filtering out samples with less than",args$normalize_count,"reads\n")
tokeep = newkey$count >= args$normalize_count
data.filt = data.norm[,tokeep]
cat("\tRemoved",sum(!tokeep),"samples\n")

# Filter OTUs based on % with zeroes
args$max_missing = 1-args$min_count_across_samples	# Late change to arguments; easiest to just put this here and not risk breaking later code
cat("Filtering out OTUs with more than",args$max_missing,"taxa with zero reads\n")
zero.fraction = apply(data.filt, MARGIN=1, FUN=function(x){sum(x==0, na.rm=T) / length(x)})
taxa.tokeep = zero.fraction <= args$max_missing
data.filt = data.filt[taxa.tokeep,]
cat("\tRemoved",sum(!taxa.tokeep, na.rm=T),"taxa\n")

# Print out results
cat("\tOriginal data file had",nrow(data),"microbial taxa across",ncol(data),"samples\n")
cat("\tFiltered data file has",nrow(data.filt),"microbial taxa across",ncol(data.filt),"samples\n")
cat("Outputting new data to",args$outfile,"\n")
data.filt = data.frame(otu=rownames(data.filt), data.filt)
write.table(data.filt, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T)

# Print out diagnostic plots
cat("Outputting grpahics to",args$printfile,"\n")
breaks=25
png(args$printfile, width=500, height=1200)
  par(mfrow = c(3,1), cex=1.2)
  # Pre-filter sample depth
  hist(newkey$count, main="Pre-filtering distribution of sample read depths", xlab="Sample depth", ylab="count", col="blue")
  legend(x="topleft", legend=paste(ncol(data),"samples"))
  
  # Pre-filter OTU/taxa zerp counts
  hist(zero.fraction, main="Pre-filtering distribution of microbial 0-fraction", xlab="Fraction 0s among samples", ylab="count", col="red")
  legend(x="topleft", legend=paste(nrow(data),"OTU taxa"))
  
  # Post-filter OTU/taxa zerp counts
  hist(zero.fraction[zero.fraction <= args$max_missing], main="Post-filtering distribution of microbial 0-fraction", xlab="Fraction 0s among samples", ylab="count", col="red")
  legend(x="topleft", legend=paste(nrow(data.filt),"OTU taxa"))
dev.off()
