#! /usr/bin/Rscript

# Do clustering as per geneSLOPE (doi: 10.1534/genetics.116.193987), but without their software because it is pretty inflexible

library(argparse)
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input file of hits by trait (from 4w_CombineParsedGwasResults.r)")
parser$add_argument("-l", "--ld", help="TASSEL linkage disequilibrium output file that includes pairwise LD for all SNP hits")
parser$add_argument("-g", "--genos", help="Hapmap file of genotypes used to make the LD file (used to match marker position to name)")
parser$add_argument("-c", "--cutoff", type="double", default=0.2, help="Cutoff for SNPs to be considered in different linkage blocks")
parser$add_argument("-o", "--outprefix")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/282MaizeLeaf16s_cleaned/4_GWAS/")
# args=parser$parse_args(c("-i","4w_combined_gwas.combined.txt", "-l", "4w_combined_gwas.hits.ld.txt", "-g", "4w_combined_gwas.hits.hmp.txt.gz", "-o","99_tmp"))

# Read in data
cat("Clustering SNPs into linkage blocks by trait for results in",args$infile,"\n")
data=read.delim(args$infile)
ld=read.delim(args$ld)
genos=read.delim(args$genos)[,1:4]	# Take only first 4 columns of metadata


# Sanity check to make sure genotype and LD match up
# # First make a key of site indices, loci, and positions
site1_key = subset(ld, select=c("Locus1", "Position1", "Site1"))
site2_key = subset(ld, select=c("Locus2", "Position2", "Site2"))
names(site1_key) = sub(names(site1_key), pattern="1$", repl="")
names(site2_key) = sub(names(site2_key), pattern="2$", repl="")
sitekey = unique(rbind(site1_key, site2_key))
sitekey = sitekey[order(sitekey$Site),]
sitekey$Site = as.numeric(sitekey$Site)	# So is the same class as the geno$index, below
# # Now make sure everything matches up
genos$index = (1:nrow(genos)) -1
if(nrow(genos) != nrow(sitekey)){stop("Different number of sites in LD file and genotype file.")}	# Check that number of sites the same
if(!identical(genos$index, sitekey$Site)){stop("Unable to match gentotypes to LD based on site index")}
if(!identical(genos$chrom, sitekey$Locus)){stop("Unable to match gentotypes to LD based on chromosome locations")}
if(!identical(genos$pos, sitekey$Position)){stop("Unable to match gentotypes to LD based on positions within chromosomes")}
# # Assuming we made it this far, make the key so that can match site indices to site names
sitekey$marker = genos$rs.


# Build matrix of LD; probably a faster way using functions, but the loop works well enough
ld_vals = matrix(NA, nrow=nrow(genos), ncol=nrow(genos), dimnames=list(genos$rs., genos$rs.))
for(i in 1:nrow(ld)){
  row=ld$Site1[i]+1	# have to add 1 because R is 1-based while TASSEL and java are 0-based
  col=ld$Site2[i]+1
  rsq=ld$R.2[i]
  if(!is.na(ld_vals[row, col]) || !is.na(ld_vals[row, col])){
	warning("Existing LD value for sites",row-1,"and",col-1,"of",ld_vals[row,col],"or",ld_vals[col,row],"will be replaced with",rsq)
  }
  ld_vals[row, col] = ld_vals[col, row] = rsq
}
# Fill in diagonal
for(i in 1:nrow(ld_vals)){ld_vals[i,i]=1}
  

# Function to do the actual clustering work
cluster_snps = function(orig, ld_vals, mincor){
  # Set up data structures
  clusters=list()
  bestsnps=list()
  targets = orig
  orig$cluster=NA
  orig$best=FALSE
  myld = subset(ld_vals, subset=rownames(ld_vals) %in% targets$Marker, select=colnames(ld_vals) %in% targets$Marker)
  
  # Actual clustering algorithm, which just calls anything with an LD greater than 'cutoff' as part of the same cluster. (Works for clusters of size 1 as well)
  i=1
  while(nrow(myld)>0 & sum(is.na(myld)) < length(myld)){	# Only stop when the LD table is empty or only full of NAs
	# Find and record best SNP remaining, based on empirical p-value and original p-value. If tied, take the first one in the genome
	tmp = targets[order(targets$empirical_pval, targets$p, targets$Chr, targets$Pos),]
	best = which(targets$Marker == tmp$Marker[1])
	
	# Record which is the best SNP
	target_snp = targets$Marker[best]
	orig$best[orig$Marker == target_snp] = TRUE
		
	# Identify all SNPs correlated with it
	tocheck = which(rownames(myld)==target_snp)
	if(length(tocheck)>1){cat("\tWARNING! Found more than one match for",target_snp,"\n")} # Sanity check
	correlated = which(abs(myld[tocheck,]) >= mincor)	# Important to use absolute value; missed that the first time.
	corr_snps = rownames(myld)[correlated]
	if(!target_snp %in% corr_snps){cat("\tWARNING! Target SNP is not correlated with itself!\n")}	# Sanity check
	
	# Save cluster
	orig$cluster[orig$Marker %in% corr_snps]=i
	
	# Remove cluster from correlations and targets
	myld=myld[-correlated, -correlated, drop=F]
	targets = subset(targets, !targets$Marker %in% names(correlated))
	
	# Increment counter
	i=i+1
  }
  return(orig)
}

# Do actual clustering by trait
mydata = subset(data, data$Marker %in% rownames(ld_vals))
mydata=split(mydata, mydata$Trait)
clusters = lapply(mydata, cluster_snps, ld_vals=ld_vals, mincor=args$cutoff)

# Report on number of clusters per trait
trait_out=paste(args$outprefix, ".clusters_per_trait.txt", sep="")
cat("\tWriting summary of clusters per trait to",trait_out,"\n")
traitsummary = sapply(clusters, function(x){max(x$cluster)})
traitsummary = data.frame(clusters=traitsummary, trait=names(traitsummary))
write.table(traitsummary, file=trait_out, sep='\t', quote=F, row.names=F, col.names=T)

# Write out the entire combined dataset
combined_out = paste(args$outprefix, ".all.txt", sep="")
cat("\tWriting full cluster assigments to",combined_out,"\n")
combined = do.call(rbind, clusters)
write.table(combined, file=combined_out, sep='\t', quote=F, row.names=F, col.names=T)

# Now just the best hits in each cluster
best_out = paste(args$outprefix, ".best.txt", sep="")
cat("\tWriting best markers per cluster to",best_out,"\n")
best = subset(combined, combined$best)
write.table(best, file=best_out, sep='\t', quote=F, row.names=F, col.names=T)
