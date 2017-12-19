#! /usr/bin/env Rscript

library(argparse)
library(plotrix) # for rescale function
library(corrplot)	# For correlation plottings
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input table of alpha diversity stats")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
parser$add_argument("-c", "--correlation-anchor", help="Which alpha diversity metric to use as the correlation anchor for flipping negatively correlated correlations")
parser$add_argument("-p", "--prune", type='double',  default=0.99, help="Prune metrics so that no pair has a (absolute value) correlation equal to or above this")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/282MaizeLeaf16s_cleaned/9_PrettyGraphics/9f_AlphaDiversity/')
# args=parser$parse_args(c("-i","9f_alpha_diversity_normalized_data.txt","-o",'99_tmp','-m', '../../2_QiimeOtus/2f_otu_table.sample_filtered.no_mitochondria_chloroplast.key.tsv', '-f', 'subpop', 'date', '-c', 'berger_parker_d'))

# Load data
cat("Comparing alpha diversity metrics for",args$infile,"\n")
alpha=read.delim(args$infile, row.names=1)

# Rescale metrics to 0-1 range
cat("\tRescaling metrics to 0-1 range\n")
metrics = names(alpha)
rescaled = alpha
for(metric in metrics){
  rescaled[[metric]] = rescale(rescaled[[metric]], newrange=c(0,1))
}

# Flip any that are negatively correlated with the anchor metric so all are basically showing the same thing
cat("\tFlipping metrics to most consistent correlations\n")
anchor = metrics[1]
if(!is.null(args$correlation_anchor)){anchor = args$correlation_anchor}
flipped = rescaled
for(metric in metrics){
  if(cor(rescaled[[metric]], rescaled[[anchor]]) < 0){
	flipped[[metric]] = 1 - flipped[[metric]]
  }
}

# Prune to just low-correlated metrics; uses a greedy algorithm, so order of inputs can matter if two metrics are perfectly correlated
cat("\tPruning metrics to ones that are less than",args$prune,"correlated\n")
pruned = cor(flipped)
removed_metrics = character()
while(sum(pruned >= args$prune) > nrow(pruned)){
  submetrics = colnames(pruned)
#   cat("\t\t\tPruned has dimensions",dim(pruned),"\n")
  for(m in submetrics){
	if(sum(pruned[,m] >= args$prune)==1){next}	# If this metric has no other metrics that are highly correlated with it, go on to the next. (always 100% correlated with itself)
	targets = submetrics[pruned[,m] >= args$prune]
	mean_correlation = sapply(targets, function(x){mean(pruned[,x])})
	tokeep = which.min(mean_correlation)
	toremove = targets[-tokeep]
	pruned = pruned[!rownames(pruned) %in% toremove, !colnames(pruned) %in% toremove]	# Remove rows/columns
	removed_metrics = c(removed_metrics, toremove)
	break	# Jump out of loop so don't cause issues with trying to reference metrics that are no longer there
  }
}
cat("\t\t",length(removed_metrics),"metrics removed:",removed_metrics,"\n")
cat("\t\t",nrow(pruned),"metrics remaining after pruning:",colnames(pruned),"\n")
pruned = subset(flipped, select = names(flipped) %in% colnames(pruned))	# Back to data frame, not correlation matrix


# Helper function for plotting correlations
plot.cors=function(x, ...){
  # Original data
  cors=cor(x)
  corrplot.mixed(cors, lower='ellipse', upper='number')
  # Clustered
#   clustered = hclust(dist(t(x)))
  mydist = as.dist(1-abs(cor(x)))	# Most reliable way of turning correlations into distances, from http://research.stowers.org/mcm/efg/R/Visualization/cor-cluster/index.htm. (already did most of this by flipping things)
  clustered = hclust(mydist)	
  new_order = clustered$order
  corrplot.mixed(cors[new_order,new_order], lower='ellipse', upper='number', ...)
  
}

# Plot correlation matrices, both in original order and hierarchically clustered
for(i in 1:2){
  if(i==1){png(paste(args$outprefix, ".correlations.png", sep=""), width=24, height=16, res=150, units='in')
  }else{svg(paste(args$outprefix, ".correlations.svg", sep=""), width=24, height=16)}
  par(mfcol=c(2,3), cex=0.5)
  # Original data
  plot.cors(alpha, main="Original Data")
#   plot.cors(rescaled, main="Rescaled")	# Irrelevant
  plot.cors(flipped, main="Flipped")
  plot.cors(pruned, main="Pruned")
  dev.off()
}

# Write out pruned metrics to keep
output = subset(alpha, select=names(alpha) %in% names(pruned))
write.table(output, file=paste(args$outprefix, ".non_redundant_metrics.txt", sep=""), sep='\t', quote=F, row.names=T, col.names=T)

