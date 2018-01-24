#! /usr/bin/env Rscript

library(argparse)
library(gplots)
parser=ArgumentParser()
parser$add_argument("-a", "--set1", help="Bacterial dataset 1")
parser$add_argument("-b", "--set2", help="Bacterial dataset 2")
parser$add_argument("-o", "--outprefix", help="Output file")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/282MaizeLeaf16s_cleaned/4_GWAS/')
#args=parser$parse_args(c("-a","4m_diversity.txt", "-b", "4m_otus.txt", "-o","../99_tmp"))

# Load data
cat("Loading data to do correlations between\n")
a = read.table(args$set1, row.names=1, header=T)
b = read.table(args$set2, row.names=1, header=T)

# Unify and make sure samples are same
mya=t(a)
myb=t(b)
samples = intersect(rownames(mya), rownames(myb))
cat("\tIdentified",length(samples),"common samples between datasets\n")
mya = mya[match(samples, rownames(mya)),]
myb = myb[match(samples, rownames(myb)),]
if(!identical(rownames(mya), rownames(myb))){
  cat("\t\tWARNING!!! Row names do not match even after attempted matching!!!")
}

# Correlate
cors = cor(mya, myb, use='pairwise')
rsq = cors * cors

# Format for output
rsq = t(rsq)
write.table(rsq, file=paste(args$outprefix, ".raw.txt", sep=""), sep='\t', row.names=T, col.names=T, quote=F)

# Make graphic for easy lookup
na_count = is.na(rsq)
heatsq = rsq[rowSums(na_count)<nrow(na_count), colSums(na_count)<ncol(na_count)]
png(paste(args$outprefix, ".raw.png", sep=""), width=10, height=20, units='in', res=300)
  heatmap.2(heatsq, trace='none', margins=c(15, 30))
dev.off()

# Make prettier/easier to see
rsq=as.data.frame(rsq)
top10 = lapply(1:ncol(rsq), function(i){
  mycors = rsq[[i]]
  names(mycors) = rownames(rsq)
  mycors = sort(mycors, decreasing=T)[1:10]
  paste(names(mycors), "( r2 = ", mycors,")")
})
top10=do.call(cbind, args=top10)
colnames(top10)=names(rsq)
write.table(top10, file=paste(args$outprefix, ".top10.txt", sep=""), sep='\t', row.names=F, col.names=T, quote=F)

# Different format of making it easier to see
long = lapply(1:ncol(rsq), function(i){
  mycors = rsq[[i]]
  names(mycors) = rownames(rsq)
  mycors = sort(mycors, decreasing=T)[1:10]
  data.frame(metric=names(rsq)[i], rank=1:length(mycors), r2=round(mycors, 3), clade=names(mycors))
})
long=do.call(rbind, long)
long$r2 = format(long$r2, nsmall = 3)	# Force to show 3 decimal places
write.table(long, file=paste(args$outprefix, ".top10_long.txt", sep=""), sep='\t', row.names=F, col.names=T, quote=F)