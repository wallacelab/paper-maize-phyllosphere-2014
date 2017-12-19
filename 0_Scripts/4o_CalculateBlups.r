#! /usr/bin/Rscript

# Convert OTU counts (transformed or otherwise) into BLUPs
# Goal of this script is to hand in raw data and get BLUPs out (1 per genotype). Meant to be modular so I can add additional steps later to check how they affect things

library(argparse)
library(lme4)	
library(parallel)
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input file of traits as a matrix")
parser$add_argument("-o", "--outfile", help="Output file of traits as a matrix")
parser$add_argument("-k", "--keyfile", help="QIIME-formatted key file of sample metadata")
parser$add_argument("-t", "--target-column", help="Column in keyfile that indicates the individual genotypes/lines to get BLUPs for")
parser$add_argument("-c", "--covariates", nargs="*", help="Name of the columns in the keyfile to use as covariates")
parser$add_argument("-n", "--num-cores", default=8, type="integer", help="Number of parallel cores to run")
parser$add_argument("-s", "--seed", default=1, type="integer", help="Random seed for selecting traits to output")
parser$add_argument("-p", "--plotfile", help="Output file to print distributions for")
parser$add_argument("-l", "--log-compare", default=FALSE, action="store_true", help="Compare models with/without log transform and take the one with better model fit")
parser$add_argument("--num-traits-to-plot", type="integer", default=10, help="Number of traits to plot distributions for")
parser$add_argument("--debug", default=FALSE, action="store_true", help="Do debug mode")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/282MaizeLeaf16s_cleaned/4_GWAS/')
# args=parser$parse_args(c("-i","4m_diversity.with_flowering.txt","-o",'99_tmp.txt','-k','../2_QiimeOtus/2f_otu_table.sample_filtered.no_mitochondria_chloroplast.key.tsv','-c',"tissue_date", "collector", '--debug', '-p', '99_tmp.png','-t','Description','--log-compare'))

cat("Using mixed-model regression to collapse count data to BLUPs\n")
data=read.delim(args$infile, check.names=F, row.names=1)
key=read.delim(args$keyfile, check.names=F)
target_column=args$target_column

# Build data frame with the appropriate setup
transdata=t(data)
transmatch = match(key$'#SampleID', rownames(transdata))
mydata = suppressWarnings(data.frame(key, transdata[transmatch,,drop=F]))
mydata=subset(mydata, mydata$X.SampleID %in% names(data))	# Trim to plots I actually have OTU data for

# convert strings to factors for requires columns
for(c in c(args$covariates, target_column)){
  if(class(mydata[,c]) == "character"){
	mydata[,c] = as.factor(mydata[,c])
  }
}

# Run the model
effect.model=paste("(1|",c(args$covariates, target_column),")", sep="", collapse=" + ")	# Covariates as fixed effects, plot as random effect
cat("\tStatistical model is ~",effect.model,"\n")
traits = rownames(data)
if(args$debug){traits=traits[1:5]}
cat("Running regression on", length(traits),"traits across",args$num_cores,"computer cores. (This may take a while.)\n")

model.fun=function(t, algorithm, log=FALSE){  
  mymod = formula(paste(t, effect.model, sep="~"))
  subdata	= subset(mydata, !is.na(mydata[,t]))# Hack b/c of missing data in flowering time
  if(log){
	subdata[,t] = suppressWarnings(log(subdata[,t]))
	infinite = !is.finite(subdata[,t])
	subdata[infinite,t] = NA
  }
  result = lmer(mymod, data=subdata)
  return(result)
}

models = mclapply(traits, mc.cores=args$num_cores, function(t){
  tryCatch(model.fun(t, algorithm=args$algorithm), warning=function(x){return(NA)}, error=function(x){return(NA)})
})
names(models)=traits

if(args$log_compare){
  log_models = mclapply(traits, mc.cores=args$num_cores, function(t){
	tryCatch(model.fun(t, algorithm=args$algorithm, log=TRUE), warning=function(x){return(NA)}, error=function(x){return(NA)})
  })
  names(log_models)=paste("log",traits,sep="_")
  
  residual_variance=function(m){
	if(class(m)=="logical"){return(99)}	# For NA results
	components = as.data.frame(VarCorr(m))
	variance = components$vcov
	resid_fract = variance[components$grp=="Residual"] / sum(variance)
	return(resid_fract)
  }
  orig_residuals = sapply(models, residual_variance)
  log_residuals = sapply(log_models, residual_variance)
  
  orig_better = orig_residuals < log_residuals  
  best_models = ifelse(orig_better, yes=models, no=log_models) # Take the model that has lower fraction residual variance
  names(best_models) = ifelse(orig_better, yes=names(models), no=names(log_models)) # Take the model that has lower fraction residual variance
  
  models = best_models
}


# Remove models that gave warnings or errors
toremove = is.na(models)
badmodels = models[toremove]
goodmodels = models[!toremove]
cat("Removed",length(badmodels),"out of", length(toremove),"traits whose models gave warnings or errors.",length(goodmodels),"traits remain to get BLUPs from.\n")

blups=mclapply(goodmodels, function(m){
  # Get BLUPs for each plot
  blups = data.frame(ranef(m)[[target_column]])
  names(blups)[1] = rnorm(1)	# Random name to prevent clashes when merging in the next step
  return(blups)
})

blups=Reduce(function(x,y){
  tmp=merge(x,y, by="row.names", all=T)
  rownames(tmp) = tmp$Row.names
  tmp$Row.names=NULL
  return(tmp)
  }, blups)
names(blups) = names(goodmodels)

# Write BLUPs in TASSEL format
write("<Phenotype>", file=args$outfile)
write(c("taxa", rep("data", ncol(blups))), file=args$outfile, append=T, ncol=ncol(blups)+1, sep='\t')
blups = data.frame(rownames(blups), blups)
names(blups)[1]="Taxon"
write.table(blups, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T, append=T)


# plot out distributions before and after
traits=names(goodmodels)	# Redfine traits so only use those that have a before-and-after
nplots = args$num_traits_to_plot
if(length(traits) < nplots) { nplots = length(traits) }	# Make sure don't go over
cat("Printing distributions for",nplots,"traits before and after BLUP creation\n")
set.seed(args$seed)
toplot = sort(sample(traits, size=nplots, replace=F))
nbreaks=20
png(args$plotfile, width=1000, height=500*nplots)
  par(mfrow=c(nplots, 2), cex=1.2)
  for(t in toplot){
	orig_name = sub(t, pattern="^log_", repl="")	# If was log transformed
	raw=as.numeric(data[rownames(data)==orig_name,])
	hist(raw, breaks=nbreaks, col="blue", main=paste("Raw data distribution for",t))
	hist(blups[,t], breaks=nbreaks, col="darkgreen", main=paste("BLUP distribution for",t))
  }
dev.off()