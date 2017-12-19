#! /usr/bin/env Rscript

library(argparse)
parser=ArgumentParser()
parser$add_argument("--orig", help="QIIME PC file from original data")
parser$add_argument("--norm", help="QIIME PC file from normalzied data")
parser$add_argument("-o", "--outfile", help="Output PNG of results")
parser$add_argument("-n", "--num_pcs", type="integer", default=10, help="Number of PCs to compare")
parser$add_argument("--name", default="unknown", help="Name for this comparison set (usually the distance metric used)")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/282MaizeLeaf16s_cleaned/9_PrettyGraphics/9h_BetaDiversity/')
#args=parser$parse_args(c("--norm","9g_qiime_diversity_normalized/bray_curtis_pc.txt","--orig", "9g_qiime_diversity_orig/bray_curtis_pc.txt", "-o","99_tmp/png", "--name", "Bray-Curtis"))

# Load data
cat("Comparing beta diversity for normalized and original data\n")
read_pcs=function(infile){
  data = read.delim(infile, skip=9, row.names=1, header=F)
  data = subset(data, !rownames(data) %in% c("Biplot", "Site constraints"))
  colnames(data) = paste("PC", 1:ncol(data), sep="")
  return(data)
}
orig=read_pcs(args$orig)
norm=read_pcs(args$norm)

# Subset down to common samples
samples = intersect(rownames(orig), rownames(norm))
orig=orig[match(samples, rownames(orig)),]
norm=norm[match(samples, rownames(norm)),]
if(!identical(rownames(orig), rownames(norm))){
  warning("Row names do not match!")
}

cors = sapply(1:args$num_pcs, function(x){
  abs(cor(orig[,x], norm[,x]))	# Since sign of PCs is somewhat arbitrary, just absolute-value it for convenience
})

# Plot
png(args$outfile)
  plot(cors, col="blue", pch=20, cex=2, xlab="PC#", ylab="Abs(correlation)", main=args$name, ylim=c(0,1))
  text(x=1:args$num_pcs, y=cors-0.1, labels=round(cors, 3))
dev.off()
