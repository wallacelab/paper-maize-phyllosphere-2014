#! /usr/bin/Rscript

library(argparse)
library(KEGGREST)
parser=ArgumentParser()
parser$add_argument("-i", "--infile")
parser$add_argument("-o", "--outfile")
parser$add_argument("-n", "--num-accessions-at-once", type="integer", default=10, choices=seq(1,10), help="Number of accessions to retrieve at once. Database limits it to at most 10")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/282MaizeLeaf16s_cleaned/4_GWAS/")
# args=parser$parse_args(c("-i","4t_kegg_accessions.txt","-o","99_tmp.txt"))

# Load data and split into groups
terms = scan(args$infile, what=character())
cat("Splitting",length(terms),"KO terms into groups of size",args$num_accessions_at_once,"\n")
groups = ceiling(seq_along(terms)/args$num_accessions_at_once)
terms = split(terms, groups)

# Retrieve KO summaries
annotations = lapply(terms, keggGet)
annotations=do.call(c, annotations)

# Make final dataframe to output
outterms = sapply(annotations, function(x){x$ENTRY})
outnames = sapply(annotations, function(x){paste(x$NAME, x$DEFINITION, sep=": ")})
outfuncs = sapply(annotations, function(x){
  f1 = paste(x$PATHWAY, collapse="; ")
  f2 = paste(x$DBLINKS, collapse="; ")
  f = paste(f1, "(synonym", f2, ")")
  return(f)
})
output = data.frame(ko=outterms, name=outnames, func=outfuncs)
write.table(output, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T)

cat(nrow(output),"annotations written to",args$outfile,"\n")