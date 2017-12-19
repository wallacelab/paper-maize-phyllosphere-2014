#! /usr/bin/Rscript

# Take my plate key and turn it into a QIIME key. Also do a basic ANOVA analysis to see if I can find variables confounded with depth

options(stringsAsFactors=F)
library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infile")
parser$add_argument("-m", "--mapfile", help="Mapping file to write out")
# parser$add_argument("-a", "--anova-out", help="Output report of ANOVA on depth by other variables")
# parser$add_argument("-r", "--readfile", help="READ-DISTRIBUTION.txt file from Minimum Entropy Composition for determining read distributions")
parser$add_argument("-l", "--linekey", help="Key of inbreds and their subpops (stiff-stalk, non-stiff-stalk, tropical, etc.)")
parser$add_argument("--raw-depth-only", default=FALSE, action="store_true", help="Whether to output only raw read depth instead of also MED depth")
args=parser$parse_args()
#setwd("/media/STORAGE/Working_Files/Metagenomics/16sLeafAmplification/20150601_Batch1_Miseq/")
# args=parser$parse_args(c("-i","0_plate_key.txt","-d","1_raw_read_counts.txt","-m","99_map","-r","99_anova", "-l", "0_inbred_subpops.txt"))

# Read in data
src=read.delim(args$infile)
# reads=read.delim(args$readfile, row.names=1)
inbreds=read.delim(args$linekey)

# Convert columns
key=data.frame(SampleID=gsub(src$Full.Label, pattern='_', repl='.'))
key$LinkerPrimerSequence = key$BarcodeSequence = ''
key$tissue_date=src$Tissue.collected
key$time = ifelse(grepl(key$tissue_date, pattern='LMAD'), yes='day', no='night')
key$rna_plate = src$Box.No.
key$dna_plate=src$DNA_Plate
key$row=src$Row..extrapolated.
key$col=src$Column..extrapolated
key$plate_row=paste(key$dna_plate, key$row, sep='.')
key$plate_col=paste(key$dna_plate, key$col, sep='.')
key$collector=src$Collector
key$date=src$Date.of.collection

# Add inbred subpopulation; requires some manipulations to get it to match the file
src$inbred=sub(src$Accession, pattern=" ", repl="")
src$inbred=toupper(src$inbred)
inbreds$Inbred = toupper(inbreds$Inbred)
key$subpop = inbreds$Subpopulation[match(src$inbred, inbreds$Inbred)]

# Add depth
# rownames(reads)=gsub(rownames(reads), pattern="_", repl=".")
# readmatch = match(key$SampleID, rownames(reads))
# key$raw_depth = rowSums(reads)[readmatch]
# if(!args$raw_depth_only){
#   key$MED_depth = reads$represented_reads[readmatch]
# }

# # # # # Add Env and SourceSink for use in SourceTracker, with blanks set to "sources" and everything else a "sink"
# # # # key$Env=key$SampleID
# # # # 	key$Env[grep(key$SampleID, pattern="BLANK")]="blank"
# # # # key$SourceSink="sink"
# # # # 	key$SourceSink[grep(key$SampleID, pattern="BLANK")]="source"
# # # # 	


# Write out new key
key$Description=gsub(src$Accession, pattern='[#/)(]', repl='')	# Has to be the last column
names(key)[1]='#SampleID'
write.table(key, file=args$mapfile, sep='\t', quote=F, row.names=F, col.names=T)

# # Do basic ANOVA analysis
# if(!is.null(args$anova_out)){
#   model = aov(raw_depth ~ time + tissue_date + collector + date + rna_plate + plate_row + plate_col, data=key)
#   write("## Sequential sum-of-squares test for model fit (= significant when added in order) for raw read depth ##\n", file=args$anova_out)
#   capture.output(anova(model), file=args$anova_out, append=T)
#   write("\n\n## Drop-one sum-of-squares test for model fit ( = significant when added last) ##\n", file=args$anova_out, append=T)
#   capture.output(drop1(model, scope=~., test="F"), file=args$anova_out, append=T)
# 
# 
#   model = aov(MED_depth ~ time + tissue_date + collector + date + rna_plate + plate_row + plate_col, data=key)
#   write("\n\n## Sequential sum-of-squares test for model fit (= significant when added in order) for MED node depth ##\n", file=args$anova_out, append=T)
#   capture.output(anova(model), file=args$anova_out, append=T)
#   write("\n\n## Drop-one sum-of-squares test for model fit ( = significant when added last) ##\n", file=args$anova_out, append=T)
#   capture.output(drop1(model, scope=~., test="F"), file=args$anova_out, append=T)
# }