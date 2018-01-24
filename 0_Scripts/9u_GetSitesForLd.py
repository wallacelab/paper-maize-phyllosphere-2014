__author__ = 'jgwall'

import argparse
import gzip
import pandas as pd

debug = False


def main():
    args = parse_args()
    print("Getting marker information for sites in",args.infile)

    # Load marker names
    hits = pd.read_table(args.infile)
    markers = pd.unique(hits['markers'].dropna())
    markers = sorted(markers)

    # Make key of marker names and indices
    i, key = 0, dict()
    for line in open(args.markers):
        line = line.strip().split('\t')
        mymarker=line[0]
        if mymarker in key: print("\tWarning! duplicate marker names:",mymarker)
        key[mymarker] = i
        i+=1

    # Get indices for each marker
    indices = [key[m] for m in markers]

    # Write out
    output = [str(m) + "\t" + str(i) + "\n" for m,i in zip(markers, indices)]
    open(args.outfile, 'w').writelines(output)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-m", "--markers", help="File of marker names (= the first few columns of the Hapmap file of genotypes that will be used for LD calculations)")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

if __name__ == '__main__': main()