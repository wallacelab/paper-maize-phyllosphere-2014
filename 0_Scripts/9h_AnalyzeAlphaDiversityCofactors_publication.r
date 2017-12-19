#! /usr/bin/env Rscript

library(argparse)
library(car)	# For Type II ANOVA
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input table of alpha diversity stats")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
parser$add_argument("-m", "--mapfile", help="QIIME mapping file")
parser$add_argument("-f", "--factors", nargs="*",  default="Description", help="Which factors to test for differences")
parser$add_argument("-n", "--num_factors", type="integer", help="Number of hidden factors to use for analysis")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/282MaizeLeaf16s_cleaned/9_PrettyGraphics/')
# args=parser$parse_args(c("-i","9f_alpha_diversity_normalized_data.txt","-o",'99_tmp','-m', '../2_QiimeOtus/2f_otu_table.sample_filtered.no_mitochondria_chloroplast.key.tsv', '-f', 'subpop', 'date', '--hidden-factors'))

# Load data
cat("Analyzing alpha diversity for",args$infile,"\n")
alpha=read.delim(args$infile, row.names=1)
key=read.delim(args$mapfile, row.names=1)

# Match up to key
keymatch = match(rownames(alpha), rownames(key))
key = key[keymatch,]
data = data.frame(alpha, key)

# Helper function to run analysis and write output
get_significance=function(data, cofactors, outprefix, title, col){
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
  
  # Plot for output for easy check
  outpng = paste(outprefix, ".pvals.png", sep="")
  ss=as.data.frame(t(ss))
  ss$Residuals=NULL
#   png(outpng, width=500, height=1000)
# 	par(mfrow=c(2,1))
	boxplot(-log10(pvals), ylab=expression(bold("-log"[10]*" p-values")), col=col, main=title)
	boxplot(ss, ylab="Faction total variance", col=col)
#   dev.off()
}



  

# Set up hidden factor analysis 
cat("Peforming hidden factor analysis (=principal components of the cofactors)\n")
cofactor_set=subset(data, select=names(data)%in%args$factors)
cofactor_set = as.data.frame(lapply(cofactor_set, as.factor))
cofactor_formula = formula(paste("dummy", paste(args$factors, collapse=" + "), sep=" ~ "))
dummy=rep(0, nrow(cofactor_set))
cofactor_matrix=model.matrix(cofactor_formula, data=cofactor_set)[,-1]	# Remove 1st column because is the intercept

# Run hidden factor analysis  
hiddens = tryCatch({factanal(cofactor_matrix[,-1], factors=args$num_factors, scores="regression")$scores}, error=function(x){return(NA)})
if(class(hiddens)=="logical" && is.na(hiddens)){	# Check if failed
  cat("\tUnable to do hidden factor analysis with",n,"hidden factors given specified covariates\n")
}
hiddens=data.frame(hiddens)
# names(hiddens) = paste("HF",args$num_factors,names(hiddens), sep="")
myfactors = names(hiddens)
newdata=data.frame(data, hiddens)


# Get graphic
# png(paste(args$outprefix,".png", sep=""), width=6, height=6, units="in", res=150)
svg(paste(args$outprefix,".svg", sep=""), width=6, height=6) 
  par(mfcol=c(2,2), las=2, mar=c(5,4,1,1), font.lab=2, font.axis=2)
  # Get directly on selected factors
  get_significance(data=data, cofactors=args$factors, outprefix=args$outprefix, title="Raw cofactors", col="cornflowerblue")
  get_significance(data=newdata, cofactors=myfactors, outprefix=paste(args$outprefix, ".hf",args$num_factors, sep=""), title="Hidden factors", col="palegreen")
dev.off()