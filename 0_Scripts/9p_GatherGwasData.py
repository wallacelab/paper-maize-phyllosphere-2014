__author__ = 'jgwall'

import argparse
import numpy as np
import pandas as pd
import re

debug = False

def main():
    args = parse_args()
    results = [pd.read_table(r) for r in args.raw]
    print("Loaded GWAS results from", len(args.raw), "input files")

    # Remove rows not associated with markers
    for i in range(len(results)):
        r = results[i]
        results[i] = r.loc[ ~np.isnan(r['Chr']),:].copy()   # Copy so don't have issues with setting stuff on slices
        # print("\t\tBefore removing NAN had",len(r),"rows of data; afterward have",len(results[i]))
        if (len(results[i]) - len(r)) > 1:
            print("\t\tWARNING! When removing <NA> chromosones, had", len(results[i]) - len(r), "rows removed (expect 0 or 1)")

    # Confirm only have a single trait
    traits = [np.array(r['Trait']) for r in results]
    traits=np.unique(np.concatenate(traits))  # collapse
    if len(traits) != 1:
        print("\tWARNING!!! Expected 1 trait but got",len(traits),":",traits,"; results may not be interpretable")

    # Get adjusted p-values if permutations supplied
    if args.perms:
        print("\tCalculating empirical p-values based on permutation p-values in",args.perms)
        results = calc_empirical_pvals(results, args.perms)
        # Filter and print out
        print("\t\tWriting out results to",args.outfile)
        output = pd.concat(results)
        output = output.sort(columns=['Chr', 'Pos'])
        output.to_csv(args.outfile, sep='\t', index=None)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--raw", nargs="*", help="List of files with GWAS results; requires columns Trait, Chr, Pos, and p")
    parser.add_argument("-o", "--outfile", help="Output file")
    parser.add_argument("-p", "--perms", help="File of permuted GWAS results for calculating empirical p-values")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def calc_empirical_pvals(results, permfile):
    perms = pd.read_table(permfile)

    # Convert to a dictionary of sorted numpy arrays for easier handling
    permvals = dict()
    for c in perms.columns:
        newname = re.sub(string=c, pattern="^chr", repl="")
        newname = int(newname)
        perm_pvals = np.array(sorted(perms[c]))  # Sort first to make p-val calc easier
        perm_pvals = perm_pvals[~np.isnan(perm_pvals)]  # Remove NAs
        permvals[newname] = perm_pvals

    # Go through and assign empirical p-values to each one
    for r in results:
        mychrom =np.unique(r['Chr'])
        if len(mychrom) != 1:
            print("WARNING!!! Trying to get p-values by chromosome but not 1 chromosome found:", mychrom)
        mychrom = mychrom[0]    # Get out of a list
        r['empirical_pval'] = np.searchsorted(a=permvals[mychrom], v = r['p'], side='left') / len(permvals[mychrom])

    return(results)


if __name__ == '__main__': main()