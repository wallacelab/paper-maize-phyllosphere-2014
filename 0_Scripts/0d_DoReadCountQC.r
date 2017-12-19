#! /usr/bin/Rscript

# Do QC on my raw read counts 
# Args: (1) Input 2-column file of file names and # reads; (2) Output png

args=commandArgs(T)
options(stringsAsFactors=F)
#setwd('/media/STORAGE/Working_Files/Metagenomics/16sLeafAmplification/20150511_InitialCheckData/')
#args=c('1_raw_read_counts.txt','99_tmp.png')

infile=args[1]
outpng=args[2]

# Combine left & right reads
data=read.delim(infile)
# data$sample=sub(data$file, pattern="_L001_.+", repl="")
# data$sample=sub(data$sample, pattern="0_RawData//", repl="")

#Color code
data$color=rep("gray", nrow(data))
data$color[grepl(data$sample, pattern="LMAD")] = "darkgreen"
data$color[grepl(data$sample, pattern="LMAN")] = "darkblue"
data$color[grepl(data$sample, pattern="BLANK")] = "red"

data=data[order(data$reads, decreasing=T), ]
png(outpng, width=2000, height=1000)
	par(mar=c(20,4,4,2), las=2)
	barplot(data$reads, names.arg=data$sample, col=data$color, main="Raw read counts")
	legend(x="topright", fill=c("darkgreen","darkblue","red","gray"), legend=c("Day samples","Night samples","Blanks","Other"), cex=1.5)
dev.off()
