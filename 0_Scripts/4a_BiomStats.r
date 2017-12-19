#! /usr/bin/env Rscript

# Plot and parse some basic stats on the distribution of OTUs in a biom file

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i","--infile",help="Observation-wise summary of a biom file ('biom summarize-table --observations')")
parser$add_argument("-o","--outtext",help="Text file to output node cutoff data to")
parser$add_argument("--outgraph",help="Graphical file to output summary data to")
parser$add_argument("-l", "--log",choices=c("x","y","xy"), default="x", help="Which axes to log-transform")
parser$add_argument("-n","--top-n", type="integer", default=50, help="Select the top N otus/samples to output to a file")
parser$add_argument("-t","--topfile",help="File to output the top N otu/sample names to to")
parser$add_argument("-r", "--raw-cutoffs", nargs="*", type="double", default=c(0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
					help="Cutoffs to apply when looking at fraction relative to the total. (That is, samples/otus that have at least this fraction of total reads individually)")
parser$add_argument("-c", "--cum-cutoffs", nargs="*", type="double", default=c(0.5, 0.7, 0.8, 0.9, 0.95),
					help="Cutoffs to apply when looking at fraction relative to the total. (That is, samples/otus that have at least this fraction of total reads individually)")
parser$add_argument("--count-cutoffs", nargs="*", type="integer", default=c(0, 100, 500, 1000, 1500, 2000, 5000, 10000),
					help="Cutoffs to apply when looking at fraction relative to the total. (That is, samples/otus that have at least this fraction of total reads individually)")
args=parser$parse_args()

# Stock values
raw_cutoffs=args$raw_cutoffs
cum_cutoffs=args$cum_cutoffs
count_cutoffs=args$count_cutoffs
set.seed(1)
cutoff_colors=sample(rainbow(max(length(raw_cutoffs), length(cum_cutoffs), length(count_cutoffs))))

# Load data and slice out just the OTU counts
summary=scan(args$infile, sep='\n', what=character())
cutoff = which(summary=="Counts/sample detail:")
counts = summary[cutoff+1:length(summary)]
otus = sub(counts, pattern=": .+", repl="")
counts=as.numeric(sub(counts, pattern=".+: ", repl=""))
names(counts)=otus
counts=counts[!is.na(counts)]

# OTUs sorted by abundance
counts = sort(counts, decreasing=T)
fract = counts / sum(counts)

# Cumulative amount
cumfract = cumsum(counts) / sum(counts, na.rm=T)


# Number with at least X raw reads
pass_rawcount = sapply(count_cutoffs, function(c){
  sum(counts >= c)
})
pass_counts = data.frame(count_cutoff=count_cutoffs, n_nodes=pass_rawcount)


# Number of nodes with at least x% of total reads
num_pass = sapply(raw_cutoffs, function(c){
  sum(fract >= c)
})
pass_cutoff = ceiling(raw_cutoffs * sum(counts))
pass_raw = data.frame(raw_percent=raw_cutoffs, min_count=pass_cutoff, n_nodes=num_pass)


# Number of nodes to reach x% of total reads
cum_pass = sapply(cum_cutoffs, function(c){
  min(which(cumfract >= c))
})
pass_cum = data.frame(cum_percent=cum_cutoffs, min_count=counts[cum_pass], n_nodes=cum_pass)



# Output text descriptions
write("# Number with at least __ individual reads", file=args$outtext)
write.table(pass_counts, file=args$outtext, sep='\t', quote=F, row.names=F, col.names=T, append=T)
write("# Number with at least __ fraction of the total reads", file=args$outtext, append=T)
write.table(pass_raw, file=args$outtext, sep='\t', quote=F, row.names=F, col.names=T, append=T)
write("\n# Number to include at least __ fraction of the cumulative total reads", file=args$outtext, append=T)
write.table(pass_cum, file=args$outtext, sep='\t', quote=F, row.names=F, col.names=T, append=T)

# Plots
png(args$outgraph, width=1800, height=600)
  par(mfrow=c(1,3), cex=1.2)
  
  # OTU cutoffs
  plot(counts, xlab="sample or otu #", ylab="read count", main="Individual fraction cutoff", log=args$log)
	abline(h=count_cutoffs, col=cutoff_colors)
	text(x=1, y=count_cutoffs, labels=count_cutoffs, pos=3, col=cutoff_colors)
  plot(counts, xlab="sample or otu #", ylab="read count", main="Individual fraction cutoff", log=args$log)
	abline(v=num_pass, col=cutoff_colors)
	text(y=1, x=num_pass, labels=raw_cutoffs, pos=4, col=cutoff_colors)
  plot(cumfract, xlab="sample or otu #", ylab="% total read count", main="Cumulative fraction cutoffs", log=args$log)
	abline(h=cum_cutoffs, col=cutoff_colors)
	text(x=1, y=cum_cutoffs, labels=cum_cutoffs, pos=3, col=cutoff_colors)
	
dev.off()

# Output the top N OTU names to a file
if(!is.null(args$top_n) && !is.null(args$topfile)){
  tops = counts[1:args$top_n]
  write(names(tops), file=args$topfile)
}