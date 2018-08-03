__author__ = 'jgwall'

import argparse

debug = False


def main():
    args = parse_args()
    print("Getting number of meaningful sites from",len(args.infiles),"input files")


    sites = set()
    i = 0
    for infile in args.infiles:
        # Open file and parse header
        IN = open(infile, "r")
        header = IN.readline().strip().split()
        markerID, pvalID = header.index("Marker"), header.index("p")

        # Go through file and identify sites that pass filter
        for line in IN:
            data = line.strip().split('\t')
            if float(data[pvalID]) <= args.max_p:
                sites.add(data[markerID])
        IN.close()

        # Print out progress
        i += 1
        if i % 100 == 0: print("\tProcessed", i, "files")
        if debug and i > 5: break   # For debugging
    print("Identified", len(sites), "sites with max pval", args.max_p, "in at least 1 trait")

    # Write results to a file
    print("Writing results to",args.outfile)
    OUT = open(args.outfile, "w")
    for site in sorted(sites):
        OUT.write(site + "\n")
    OUT.close

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*", help="Input files of p-values from original GWAS")
    parser.add_argument("-o", "--outfile", help="Output file")
    parser.add_argument("-p", "--max-p", type=float, default=0.001, help="Maximum p-value in any trait for a site to be considered meaningful")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


if __name__ == '__main__': main()