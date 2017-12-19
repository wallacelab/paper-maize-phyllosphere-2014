#! /usr/bin/python3 

# Fix my sample names in the .fna file to use periods instead of underscores for the same name, and also add a delimiter for the read number

import argparse
import gzip

# Command-line arguments
parser=argparse.ArgumentParser()
parser.add_argument("-i","--infiles", nargs="*")
parser.add_argument("-o","--outfile")
parser.add_argument("--debug", default=False, action="store_true")
args=parser.parse_args()

print("Changing sample name underscores to periods in",len(args.infiles),"input files")

OUT = open(args.outfile, "w")
total=0
for infile in args.infiles:
    IN = gzip.open(infile, "rt") if infile.endswith(".gz") else open(infile, "r")
    n=0
    reads=dict()
    for line in IN:
        if line.startswith(">"):	# Only worry about sample names, not the DNA sequences
            n+=1
            if args.debug and n > 1000: break   # for debugging
            tokens=line.split('|')
            
            # Turn underscores to periods in names
            name=tokens[0].split("_")
            newname = ".".join(name)
            
            #Add read number to end of sample name
            if newname not in reads:
                reads[newname]=0
            readcount = reads[newname]
            reads[newname] +=1
            newname += "_" + str(readcount)
            
            # Put back into sample name
            tokens[0]=newname
            line=" ".join(tokens)
        OUT.write(line)
    IN.close()
    total+=n
    print("\t",n,"sequences from",infile,"(",total,"total )")
OUT.close()
print("Total",total,"sequences written to",args.outfile)