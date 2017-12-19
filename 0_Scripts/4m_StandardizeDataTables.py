__author__ = 'jgwall'

import argparse
import numpy as np
import re

debug = False
bad_chars = ['[', ']', ' ', '-']    # List of bad characters that will be replaced with underscores in trait names

def main():
    args = parse_args()
    print("Standardizing file",args.infile)

    # Open file handles and copy header over
    IN = open(args.infile, "r")
    OUT = open(args.outfile, "w")
    OUT.write(IN.readline())
    
    # Write original data out (unless told to skip it)
    if not args.skip_orig_data:
        # Fix trait names if requested, otherwise just copy over
        if args.fix_traitnames:
            print("Standardizing trait names to remove illegal characters")
            for line in IN:
                data=line.split('\t')
                trait=data[0]
                for c in bad_chars: # Loop over bad characters and replace with underscores
                    trait=trait.replace(c, '_')
                data[0] = trait
                OUT.write("\t".join(data))
        else:
            OUT.writelines(IN.readlines())   # Write out existing data
    IN.close()
    
    # Add flowering data if requested
    if(args.flowerfile):
        print("Adding flowering time")
        flowering = get_flowering(args.infile, args.flowerfile, args.jitter, args.jitter_deviation, args.seed)
        OUT.write("\t".join(flowering) + "\n")  
    
    
    OUT.close()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-f", "--flowerfile", help="Two-column file with line names and flowering times")
    parser.add_argument("-j", "--jitter", default=False, action="store_true", help="Whether to jitter the flowering times so repeat samples from the same row don't have the exact same amount")
    parser.add_argument("-d", "--jitter-deviation", type=float, default=0.001, help="Standard deviation to use when jittering flowering time")
    parser.add_argument("-s", "--seed", type=int, default=1, help="Random seed for jittering flowering time")
    parser.add_argument("-n", "--fix-traitnames", default=False, action="store_true", help="Whether to fix names of traits to remove illegal characters (brackets, spaces, etc)")
    parser.add_argument("-x", "--skip-orig-data", default=False, action="store_true", help="Whether to skip the traits originally in the file. (Used to get a file of just flowering time)")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

# Cut out just the row identifier for each sample
def get_row(sample):
    mymatch = re.search("14A....", string=sample)
    if mymatch is not None:
        # print("Sample",sample,"is in row",mymatch.group(0))
        return mymatch.group(0)

def get_flowering(infile, flowerfile, jitter, jitter_deviation, seed):
    samples = open(infile).readline().strip().split('\t')
    if not samples[0].startswith("LMA"): # Remove the first column header if it's included
        samples=samples[1:]

    # Load keys
    flowerkey = dict()
    for line in open(flowerfile):
        name, time = line.strip().split('\t')
        if time == "NA": time=np.nan
        flowerkey[name]=float(time)

    # Create flowering time data
    flowering = list()
    for sample in samples:
        if get_row(sample) in flowerkey:
            flowering.append(flowerkey[get_row(sample)])
        else:
            print("Unable to find",sample)
            flowering.append(np.nan)
    flowering = np.array(flowering)

    # # Add random noise to flowering time data if requested
    if jitter:
        np.random.seed(seed)
        flowering = flowering + np.random.normal(loc=0, scale=jitter_deviation, size=len(flowering))

    # Convert to string
    flowering = ["flowering_time"] +  [str(f) for f in flowering]
    
    # Change "nan" to "NA" to be consistent with other phenotypes
    for i in range(len(flowering)):
        if flowering[i] == "nan" : 
            #print("Changing flowering at index",i)
            flowering[i] = "NA"
    
    return(flowering)

if __name__ == '__main__': main()