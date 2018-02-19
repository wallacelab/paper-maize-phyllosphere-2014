#! /usr/bin/Rscript

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infiles", nargs="*", help="List of input FASTQ files")
parser$add_argument("-o", "--outfile")
args=parser$parse_args()
#setwd("/home/jgwall/Projects/282MaizeLeaf16s_cleaned/0_RawData/")
#args=parser$parse_args(c("-i","5504_8698_23201_HJC7WADXX_LMAN_8_14A0091_TAAGGCGA_CTCTCTAT_R1.fastq.gz", "5504_8698_23201_HJC7WADXX_LMAN_8_14A0091_TAAGGCGA_CTCTCTAT_R2.fastq.gz", 
# "5504_8698_23202_HJC7WADXX_LMAN_8_14A0095_TAAGGCGA_TATCCTCT_R1.fastq.gz", "5504_8698_23202_HJC7WADXX_LMAN_8_14A0095_TAAGGCGA_TATCCTCT_R2.fastq.gz",  "-o","../99_tmp.txt"))

cat("Collecting",length(args$infiles),"file names for submission to the SRA\n")

# Get sample names
data=data.frame(sample="", files=args$infiles)
data$sample=sub(data$files, pattern=".+(LMA[ND]_[^_]+_14A.+)_[^_]+_[^_]+_R..fastq.gz", repl="\\1")	# Very complex regex to get sample name from file name
data$sample=gsub(data$sample, pattern="_", repl=".")	# To match QIIME keyfile

# Split and combine
samples = split(data, data$sample)
samples=lapply(samples, function(x){
  forward = x$files[grep(x$files, pattern="_R1.fastq.gz")]
  reverse = x$files[grep(x$files, pattern="_R2.fastq.gz")]
  data.frame(sample=unique(x$sample), forward=forward, reverse=reverse)
})

# Sanity check
if(length(table(sapply(samples, nrow))) != 1){
  cat("WARNING!! Unequal lengths found for files! Some samples with >1 Forward or Reverse sequence file\n")
}

# Combine and write out
output=do.call(rbind, samples)
cat("\tCollapsed to forward and reverese reads for",nrow(output),"samples\n")
write.table(output, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T)
