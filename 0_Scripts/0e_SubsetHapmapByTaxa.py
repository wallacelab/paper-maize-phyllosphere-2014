__author__ = 'jgw87'
"""
Cut a hapmap down to just the taxa in a specified file
"""

import argparse
import gzip
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import re

taxastart=11

def main():
    args = parse_args()
    IN = get_handle(args.infile, "rt")
    OUT = get_handle(args.outfile, "wt")

    # Load taxa list
    taxa = open(args.taxafile).readlines()
    taxa  = set([t.strip() for t in taxa])

    # Get header
    header = np.array(IN.readline().strip().split('\t'))
    tokeep = np.array([False] * len(header))
    tokeep[:taxastart] = True

    # Find taxa to keep
    keepers=0
    for i in range(len(header)):
        if header[i] in taxa:   # Since metadata already set to true, this won't affect it
            tokeep[i] = True
            keepers+=1
    print("Found",keepers,"taxa (out of",len(taxa),") to keep; processing genotypes")

    # Write out new file
    OUT.write("\t".join(header[tokeep]) + "\n")
    n=0
    for line in IN:
        n+=1
        if n % 25000 == 0: print("\tProcessed",n,"lines")
        if args.debug and n > 1000: break  # For debugging
        newdata = np.array(line.strip().split('\t'))
        OUT.write("\t".join(newdata[tokeep]) + "\n")


    IN.close()
    OUT.close()

def get_handle(file, mode):
    if file.endswith(".gz"): return gzip.open(file, mode)
    else: return open(file, mode)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-t", "--taxafile")
    parser.add_argument("-d", "--debug", default=False, action="store_true")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()