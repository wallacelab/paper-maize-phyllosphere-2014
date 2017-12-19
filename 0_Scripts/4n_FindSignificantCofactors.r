#! /usr/bin/Rscript

# Convert OTU counts (transformed or otherwise) into BLUPs
# Goal of this script is to hand in raw data and get BLUPs out (1 per genotype). Meant to be modular so I can add additional steps later to check how they affect things

library(argparse)
library(lme4)
library(parallel)
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input file of traits as a matrix")
parser$add_argument("-o", "--outfile", help="Output file of p-values for significance")
parser$add_argument("-k", "--keyfile", help="QIIME-formatted key file of sample metadata")
parser$add_argument("-c", "--covariates", nargs="*", help="Name of the columns in the keyfile to use as covariates")
parser$add_argument("-n", "--num-cores", default=8, type="integer", help="Number of parallel cores to run")
parser$add_argument("-s", "--seed", default=1, type="integer", help="Random seed for selecting traits to output")
parser$add_argument("-p", "--plotfile", help="Output file to print distributions for")
parser$add_argument("-q", "--cofactors-to-test", nargs="*", help="Output file to print distributions for")
parser$add_argument("--debug", default=FALSE, action="store_true", help="Do debug mode")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/282MaizeLeaf16s_cleaned/4_GWAS/')
# args=parser$parse_args(c("-i","4m_diversity.with_flowering.txt","-o",'99_tmp.txt','-k','../2_QiimeOtus/2f_otu_table.sample_filtered.no_mitochondria_chloroplast.key.tsv','-q',"tissue_date", "dna_plate","collector", '--debug', '-p', '99_tmp.png'))

cat("Identifying significance of potential cofactors to include\n")
data=read.delim(args$infile, check.names=F, row.names=1)
key=read.delim(args$keyfile, check.names=F)

# Pull out the plot from the sample name
key$plot = sub(key$'#SampleID', pattern=".+\\.", repl="")

# Build data frame with the appropriate setup
# rownames(data) = paste("trait",rownames(data),sep="_") # So formula cause problems with any numerically named traits
transdata=t(data)
transmatch = match(key$'#SampleID', rownames(transdata))
mydata = suppressWarnings(data.frame(key, transdata[transmatch,]))
mydata=subset(mydata, mydata$X.SampleID %in% names(data))	# Trim to plots I actually have OTU data for

# convert strings to factors for required columns
for(c in c(args$covariates, args$cofactors_to_test, "plot")){
  if(class(mydata[,c]) == "character"){
	mydata[,c] = as.factor(mydata[,c])
  }
}


# Get list of traits
traits = rownames(data)
if(args$debug){traits=traits[1:5]}

# function to run two models and compare
test.model = function(t){
	# Model specifications
	mynull = formula(paste(t, null.model, sep="~"))
	myrand = formula(paste(t, rand.model, sep="~"))
	
	# Run models
	nullresult = lmer(mynull, data=mydata, REML=FALSE)
	randresult = lmer(myrand, data=mydata, REML=FALSE)
	  
	# ANOVA
	a = anova(nullresult, randresult)
	return(a)
}

pvals=list()
for(query in args$cofactors_to_test){
  if(query %in% args$covariates){
	cat("WARNING!! Told to test", query,"but is one of the supplied covariates. Will be ignored.\n")
	next
  }
  
  cat("Testing significance of", query,"in", length(traits),"traits across",args$num_cores,"computer cores. (This may take a while.)\n")
  null.model=paste("(1|", c(args$covariates, "plot"), ")", sep="", collapse=" + ")
  rand.model=paste("(1|", c(args$covariates, query, "plot"), ")", sep="", collapse=" + ")
  cat("\tNull model is ~",null.model,"\n")
  cat("\tFull model is ~",rand.model,"\n")

  
  
  # Actual models
  anovas = mclapply(traits, mc.cores=args$num_cores, FUN=function(t){
	  tryCatch(test.model(t), warning=function(x){return(NA)}, error=function(x){return(NA)})
	})
  
  # extract pvalues
  pval_col="Pr(>Chisq)"
  
  pvals[[query]] = sapply(anovas, function(a){
	if("data.frame" %in% class(a)){
	  return(a[,pval_col][-1])
	}else{return(NA)}
  })
}
pvals=as.data.frame(pvals)

# Write out results
write.table(pvals, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T)

# Output summaries to stdout
means=sapply(pvals, mean, na.rm=T)
counts=sapply(pvals, function(x){sum(!is.na(x))})
cat("Summary of results:\n")
for(i in 1:ncol(pvals)){
  percent=round(counts[i]/nrow(pvals) * 100, digits=1)
  cat("\t",names(pvals)[i],"has mean significance",means[i],"from",counts[i],"out of",nrow(pvals),"traits (",percent,"%)\n")
}

# Plot results (log-transformed)
toplot=pvals
toplot[toplot==0] = min(toplot, na.rm=T)*0.1	# Adjust for any "pvalue = 0" issues
toplot=-log10(toplot)	
png(args$plotfile, width=800, height=800)
  par(cex=1.2)
  box = boxplot(toplot, main="Distribution of -log10 pvals")
  legend(x="topright", title="Sample Sizes", legend=paste(box$n, box$names))
dev.off()