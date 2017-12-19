__author__ = 'jgwall'

import argparse
import numpy as np
import pandas as pd
import re

debug = False


def main():
    args = parse_args()
    print("Making Annotation Enrichment signature file from",args.infile)
    data = pd.read_table(args.infile)

    # Subset to just the traits that pass empirical p-value
    good_pval = data['empirical_pval'] <= args.pval_cutoff
    good_herit = data['h2'] <= args.max_herit
    print("\tFiltering out", np.sum(~good_pval),"traits with p-values above",args.pval_cutoff)
    print("\tFiltering out", np.sum(~good_herit), "traits with heritabilities above", args.max_herit,"because are probably artifacts")
    data = data.loc[good_pval & good_herit,:]
    print("\tData filtered from",len(good_pval),"traits to",len(data),"traits")

    # Reformat traits for signature file
    traits = data['trait']
    traits = [re.sub(pattern="^log_", string=t, repl="") for t in traits]
    traits = [re.sub(pattern="_[0-9]+$", string=t, repl="") for t in traits]

    # Subset out based on set
    pattern=""
    if args.set=="ko": pattern="K"
    elif args.set=="cog": pattern="COG"
    else: print("Warning! Unable to recognize set",args.set,"to determine which temrs to filter for")
    target_traits = set()
    for t in traits:
        if t.startswith(pattern): target_traits.add(t)
    print("\t",len(target_traits),"traits identified in set",args.set)

    # Make output file
    print("Writing signature file to", args.outfile)
    OUT = open(args.outfile, "w")
    terms = ",".join(sorted(target_traits))
    OUT.write(args.name + "\t" + terms)
    OUT.close()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Parsed heritabilities from 5_AssembleMetagenomeHeritData.py")
    parser.add_argument("-s", "--set", choices=['cog','ko'], help="Whether this is COG or KEGG terms ")
    parser.add_argument("-o", "--outfile", help="Output signature file for Annotation Enrichment Analysis")
    parser.add_argument("-p", "--pval-cutoff", type=float, default=1, help="P-value cutoff for including in the signature")
    parser.add_argument("-m", "--max-herit", type=float, default=1, help="Exclude heritability above this (b/c is probably a artifact)")
    parser.add_argument("-n", "--name", help="Name to give this signature")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


if __name__ == '__main__': main()