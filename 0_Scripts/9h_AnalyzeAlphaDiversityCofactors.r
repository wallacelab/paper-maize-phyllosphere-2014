#! /usr/bin/env Rscript

library(argparse)
library(car)	# For Type II ANOVA
library(colorspace) # For plot of variance explained
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input table of alpha diversity stats")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
parser$add_argument("-m", "--mapfile", help="QIIME mapping file")
parser$add_argument("-f", "--factors", nargs="*",  default="Description", help="Which factors to test for differences")
parser$add_argument("--hidden-factors", default=FALSE, action="store_true", help="Flag to perform hidden factor analysis on the chosen cofactors (= fit their principal components)")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/282MaizeLeaf16s_cleaned/9_PrettyGraphics/9f_AlphaDiversity/')
# args=parser$parse_args(c("-i","9f_alpha_diversity_normalized_data.txt","-o",'99_tmp','-m', '../../2_QiimeOtus/2f_otu_table.sample_filtered.no_mitochondria_chloroplast.key.tsv', '-f', 'subpop', 'date', '--hidden-factors'))

# Load data
cat("Analyzing alpha diversity for",args$infile,"\n")
alpha=read.delim(args$infile, row.names=1)
key=read.delim(args$mapfile, row.names=1)

# Match up to key
keymatch = match(rownames(alpha), rownames(key))
key = key[keymatch,]
data = data.frame(alpha, key)

# Helper function to run analysis and write output
get_significance=function(data, cofactors, outprefix){
  # Go through and test factors for differences
  cat("\tCalculating Type II SS significance for",cofactors,"\n")
  cofactors = paste(cofactors, collapse=" + ")
  anovas = lapply(names(alpha), function(metric){
	  myformula = paste(metric,"~",cofactors)
	  mymodel = lm(myformula, data=data)
	  myanova = Anova(mymodel, type=2)
  })

  # Extract p-values and turn into useful data frame
  pvals = lapply(anovas, function(x){
	x = as.data.frame(x)
	x = subset(x, select=names(x) == 'Pr(>F)')
  })
  pvals = do.call(cbind, pvals)

  ss = lapply(anovas, function(x){
	x = as.data.frame(x)
	x = subset(x, select=names(x) == 'Sum Sq')
	x = x / sum(x)
  })
  ss = do.call(cbind, ss)

  # Format for output
  names(pvals) = names(alpha)
  pvals=as.data.frame(t(pvals))
  pvals$Residuals=NULL
  outfile = paste(outprefix, ".pvals.txt", sep="")
  write.table(pvals, file=outfile, sep='\t', row.names=T, col.names=T, quote=F)
  
  # Get fraction total variance explained
  fractional = 1- ss["Residuals",]
  names(fractional) = names(alpha)
  
  # Plot for output for easy check
  outpng = paste(outprefix, ".pvals.png", sep="")
  ss=as.data.frame(t(ss))
  ss$Residuals=NULL
  png(outpng, width=500, height=1000)
	par(mfrow=c(2,1))
	boxplot(-log10(pvals), xlab="cofactor", ylab="-log10 pvalues")
	boxplot(ss, xlab="cofactor", ylab="Faction total variance explained")
  dev.off()
  
  return(fractional)
}
  
# Get directly on selected factors, and save the proportion total variance explained
cofactor_fractional = get_significance(data=data, cofactors=args$factors, outprefix=args$outprefix)
rownames(cofactor_fractional)[1] = "raw_cofactors"

# Perform hidden factor analysis if requested
fractionals = list()	# List to save proportion variance explained
if(args$hidden_factors){
  cat("Peforming hidden factor analysis (=principal components of the cofactors)\n")
  cofactor_set=subset(data, select=names(data)%in%args$factors)
  cofactor_set = as.data.frame(lapply(cofactor_set, as.factor))
  cofactor_formula = formula(paste("dummy", paste(args$factors, collapse=" + "), sep=" ~ "))
  dummy=rep(0, nrow(cofactor_set))
  cofactor_matrix=model.matrix(cofactor_formula, data=cofactor_set)[,-1]	# Remove 1st column because is the intercept
  
  for(n in 1:ncol(cofactor_matrix)){
	hiddens = tryCatch({factanal(cofactor_matrix[,-1], factors=n, scores="regression")$scores}, error=function(x){return(NA)})
	if(class(hiddens)=="logical" && is.na(hiddens)){	# Check if failed
	  cat("\tUnable to do hidden factor analysis with",n,"hidden factors given specified covariates\n")
	  next
	}
	hiddens=data.frame(hiddens)
	names(hiddens) = paste("HF",n,names(hiddens), sep="")
	myfactors = names(hiddens)
	newdata=data.frame(data, hiddens)
	fractionals[[n]] = get_significance(data=newdata, cofactors=myfactors, outprefix=paste(args$outprefix, ".hf",n, sep=""))
	rownames(fractionals[[n]]) = paste("hidden_factors",n, sep="")
  }
}

# Collate proportion variance explained
variance_explained = do.call(rbind, fractionals)
if(is.null(variance_explained)){  variance_explained = cofactor_fractional  # Conditional to handle if hidden factor analysis was skipped
} else {variance_explained = rbind(cofactor_fractional, variance_explained) }

# Plot proportion variance explained
for( i in 1:2){
  if(i==1){png(paste(args$outprefix, ".variance_explained.png", sep=""), width=16, height=6, res=150, units='in')
  }else{svg(paste(args$outprefix, ".variance_explained.svg", sep=""), width=16, height=6)}
  par(mar=c(12,4,4,1), las=2, font.lab=2, font.axis=2)
  colors = rainbow_hcl(nrow(variance_explained))
  barplot(as.matrix(variance_explained), ylab="% variance explained", beside=T, legend=T, col=colors, args.legend=list(x='topleft', cex=0.6))
  dev.off()
}