#! /usr/bin/env Rscript

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infile")
parser$add_argument("-t", "--traits", nargs="*")
parser$add_argument("-o", "--outprefix")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/282MaizeLeaf16s_cleaned/6_PowerSimulation/")
# args=parser$parse_args(c("-i","../4_GWAS/4u_combined_good_traits.update_flowering.tassel.txt","-o","99_tmp", "-t", "flowering_time", "weighted_unifrac_PC1"))

# Load data
cat("Filtering phenotypes in TASSEL-formatted file", args$infile,"\n")
data = read.delim(args$infile, skip=1, header=F)
cat("\t", ncol(data)-1,"phenotypes loaded\n")

# Filter
traits = unlist(data[2,])   # Trait names are in second row of data table
targets = c("Taxa", args$traits)
tokeep = traits %in% targets
cat("\tFound",sum(tokeep),"traits to keep.\n")

# Sanity check that found all specified traits
if(sum(tokeep) != length(args$traits)+1){
    notfound = args$traits[!args$traits %in% traits]
    cat("\t\tWarning! Could not find traits",notfound,"\n")
}else{
    cat("\t\tAll traits found.\n")
}

# Subset
subset = data[,tokeep]

# Output new subsetted text file
cat("Outputting new text files\n")
for(i in 2:ncol(subset)){
    myname = subset[2,i]    # Trait name
    outfile = paste(args$outprefix, myname, "tassel.txt", sep=".")
    write("<Phenotype>", outfile)
    write.table(subset[,c(1,i)], file=outfile, row.names=F, col.names=F, quote=F, sep='\t', append=T)   # Write for just this phenotype
}

# Redo data frame more suitable for plotting histograms
subdata = data.frame(subset[-c(1:2), -1], row.names=subset[-c(1:2),1])
names(subdata)=unlist(subset[2,-1])


# Plot histograms
cat("Plotting histogram of traits\n")
png(paste(args$outprefix, ".hists.png", sep=""), width=500, height=500*ncol(subdata))
    par(mfrow=c(ncol(subdata), 1))
    tmp = lapply(names(subdata), function(trait){
        vals = as.numeric(as.character(subdata[,trait]))
        hist(vals, main=trait, breaks=30, col='darkblue', xlab="Value", ylab="Frequency", cex=1.4)
    })
dev.off()
