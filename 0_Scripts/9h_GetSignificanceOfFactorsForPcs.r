#! /usr/bin/env Rscript

library(argparse)
library(Hotelling)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="QIIME PC file")
parser$add_argument("-o", "--outprefix", help="Output file")
parser$add_argument("-k", "--keyfile", help="QIIME_keyfile")
parser$add_argument("-n", "--num_pcs", type="integer", default=3, help="Number of PCs to compare")
parser$add_argument("--include-batches", default=FALSE, action="store_true", help="Whether to compare batches based on dna_plate")
parser$add_argument("--name", default="unknown", help="Name for this comparison set (usually the distance metric used)")
parser$add_argument("--skip", default="", nargs="*", help="Any categories that should be skipped")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/282MaizeLeaf16s_cleaned/9_PrettyGraphics/9h_BetaDiversity/')
#args=parser$parse_args(c("-i","9g_qiime_diversity_normalized/bray_curtis_pc.txt", "-o","99_tmp", "--name", "Bray-Curtis", "--include-batches", "-k", "../../2_QiimeOtus/2f_otu_table.sample_filtered.no_mitochondria_chloroplast.key.tsv"))

# Load data
cat("Performing MANOVA on PCs in",args$infile,"\n")
key=read.delim(args$keyfile, row.names=1)
data = read.delim(args$infile, skip=9, row.names=1, header=F)
data = subset(data, !rownames(data) %in% c("Biplot", "Site constraints"))
colnames(data) = paste("PC", 1:ncol(data), sep="")

# Reduce down to just the top PCs
if(args$num_pcs < ncol(data)){
  cat("\tReducing to just first",args$num_pcs,"PCs\n")
  data=data[,1:args$num_pcs]
}

if(args$include_batches){
  key$batch = factor(ifelse(key$dna_plate %in% c("KAK_7_8", "KAK_9", "KAK_11", "KAK_12"), yes="batch1", no="batch2"))
}

# Match up key
keymatch = match(rownames(data), rownames(key))
newkey = key[keymatch,]
newkey=droplevels(newkey)

# Perform MANOVA
pcs=as.matrix(data)
outputs=list()
cat("Calculating significance\n")
for(label in names(newkey)){
  if(length(levels(newkey[,label])) <2){ next }
  if(label %in% args$skip) {next}	# Skip any specified by the user
  cat("\tCategory:",label,"\n")
  myform = formula(paste("pcs ~", label))
  model = manova(myform, data=newkey)
  
  # Make pairwise tests
  my_levels = levels(newkey[,label])
  combos = combn(1:length(my_levels), m=2, simplify=FALSE)
  pairwise_pvals = sapply(combos, function(c){
	test = hotelling.test(myform, data=newkey, pair=c)
	return(test$pval)
  })
  
  # Format output
  group_name=paste("###",label,"###",sep="")
  group_stats=as.data.frame(summary(model)$stats)
  group_pval = group_stats$"Pr(>F)"[rownames(group_stats)==label]
  my_output = data.frame(test=group_name, pval=group_pval)
  
  # Add in individual pairwise tests in a pretty format.
  comparisons = sapply(combos, function(x){
	  paste(my_levels[x], collapse="|")
  })
  comp_output = data.frame(test = comparisons, pval=pairwise_pvals)
  comp_output = comp_output[order(comp_output$pval),]
  my_output = rbind(my_output, comp_output)
  
  outputs[[label]] = my_output
}
  
outfile=paste(args$outprefix, ".pvals.txt", sep="")
cat("Writing output to",outfile,"\n")
output = do.call(rbind, outputs)
names(output)[2] = paste("pval",args$name,sep="_")
write.table(output, file=outfile, sep='\t', quote=F, row.names=F, col.names=T)

# Write a shorter version without the pairwise comparisons
outfile=paste(args$outprefix, ".pvals_short.txt", sep="")
suboutput = subset(output, grepl(output$test, pattern="^##"))
suboutput$test = gsub(suboutput$test, pattern="#", repl="")
write.table(suboutput, file=outfile, sep='\t', quote=F, row.names=F, col.names=T)



# Stepwise to identify effect of including other variables
potentials = names(outputs)
included = character()
cat("Performing stepwise MANOVA with",potentials,"\n")
while(length(potentials)>0){
  dependents = paste(included, collapse="+")
  if(length(dependents)==0){dependents = "1"}	# Take care of first go through
  myform = paste("pcs ~", dependents)
  
  # Calculate pvalues for adding each one in
  pvals = sapply(potentials, function(p){
	 testformula = formula(paste(myform, p, sep="+"))
	 model = manova(testformula, data=newkey)
	 stats=as.data.frame(summary(model)$stats)
	 pval = stats$"Pr(>F)"[rownames(stats)==p]
	 return(pval)
  })
  
  # Adjust model
  target = which.min(pvals)
  included = c(included, potentials[target])
  potentials = potentials[-target]
}
  
# And now that have everything in order, do simple Type I SS
dependents = paste(included, collapse="+")
myform = formula(paste("pcs ~", dependents))
model = manova(myform, data=newkey)
stats=as.data.frame(summary(model)$stats)

outfile=paste(args$outprefix, ".stepwise_manova.txt", sep="")
cat("Writing stepwise MANOVA output to",outfile,"\n")
write.table(stats, file=outfile, sep='\t', quote=F, row.names=T, col.names=T)