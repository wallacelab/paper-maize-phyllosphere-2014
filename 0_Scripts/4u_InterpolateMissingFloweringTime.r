#! /usr/bin/Rscript

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="TASSEL-formatted covariate file of flowering time")
parser$add_argument("-f", "--flowering-averages", help="File of flowering time averages for the Goodman panel from Panzea.org")
parser$add_argument("-t", "--taxa", help="File listing the taxa to try to interpolate/include")
parser$add_argument("-o", "--outfile", help="Output covariate file")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/282MaizeLeaf16s_cleaned/4_GWAS/")
# args=parser$parse_args(c("-i","4p_flowering_covariate.tassel.txt","-f", "../0_Scripts/0_goodman_average_flowering.panzea.txt", "-o","99_tmp.txt", "-t", "4t_combined_good_traits.taxa.txt" ))

# Load data
cat("Loading data to interpolate missing flowering times\n")
header=scan(args$infile, what=character(), sep='\n', nlines=2)
data=read.delim(args$infile, skip=2)
avg=read.delim(args$flowering_averages)
taxa=scan(args$taxa, what=character())

# Convert taxa to a data frame to match up
taxa = data.frame(name=taxa)

# function to standardize names
standardize=function(x){
  x=toupper(x)	# Make all uppercase
  x=sub(x, pattern=":.+", repl="")	# Remove anything after the first colon
  x=gsub(x, pattern="[-_ ]", repl="")	# Remove hyphens, underscores, and spaces
  
  # Manual hacks for the things the above didn't take care of
  x[x=="B2GOOD"]="B2"
  x[x=="W22RRSTDCS29091"] = "W22RRSTD"
  return(x)
}

data$key = standardize(data$Taxon)
avg$key=standardize(avg$taxa)
taxa$key = standardize(taxa$name)

# Test for concordance
find_missing=function(x, y){
  notfound = x[! x %in% y]
  cat("\t",length(notfound),"taxa not found:",notfound,"\n")
  cat("\t(These should be mostly CMLs)\n")
}
cat("Compare master taxa to covariate file\n")
find_missing(taxa$key, data$key)
cat("Compare covariate file to master taxa\n")
find_missing(data$key, taxa$key)
cat("Compare master taxa to file of average values\n")
find_missing(taxa$key, avg$key)

# Match up
datamatch = match(taxa$key, data$key)
avgmatch = match(taxa$key, avg$key)
combined = cbind(taxa, data[datamatch,], avg[avgmatch,])

# Interpolate via linear regression
model = lm(flowering_time ~ average, data=combined)
combined$predictions = predict(model, newdata=data.frame(taxa=combined$name, average=combined$average))
is_missing = is.na(combined$flowering_time)
combined$interp = ifelse(is_missing, yes=combined$predictions, no=combined$flowering_time)

# Format and write out
output=data.frame(Taxon=combined$name)
output[,2] = combined$interp
names(output)[2] = names(data)[2]	# Copy phenotype name over
write(header, file=args$outfile)
write.table(output, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T, append=T)


