__author__ = 'jgwall'

import argparse
import numpy as np
import pandas as pd

debug = False


def main():
    args = parse_args()
    print("Getting TASSEL site file data for traits in",args.infile,"at p-value cutoff of",args.pval_cutoff)
    data=pd.read_table(args.infile)
    sitefile = pd.read_table(args.sitefile)

    # Subset based on empirical p-value
    data=data.loc[data['empirical_pval'] <= args.pval_cutoff]

    # Make a key of trait-marker combinations
    targetkey = {str(trait) + ":" + str(marker) for trait,marker in zip(data['Trait'], data['Marker'])}
    sitekey = [str(trait) + ":" + str(marker) for trait,marker in zip(sitefile['Trait'], sitefile['Marker'])]
    targets = [hit in targetkey for hit in sitekey]

    # Sanity checks
    print("\tFound",len(data),"hits that passed p-value cutoff, which collapse to",len(targetkey),"unique trait-marker combinations (should be the same number)")
    print("\tLocated",np.sum(targets),"hits to include in output data (should also be the same)")
    siteset = set(sitekey)
    missing = [target not in siteset for target in targetkey]
    print("\t\t",np.sum(missing),"could not be found (should be 0)")

    # Print out those that are missing
    if np.any(missing):
        mymissing = list()
        for target in targetkey:
            if target not in siteset: print("\t\t\t",target)

    # Subset
    output = sitefile.loc[targets, :]

    # Output
    output = output.sort(columns=['Trait','Chr', 'Pos'])
    output.to_csv(args.outfile, sep='\t', index=False)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="File of clustered SNP hits from 4x_ClusterSnpHits.r (with the .best.txt suffix)")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-s", "--sitefile", help="TASSEL sitefile of traits and hits")
    parser.add_argument("-p", "--pval-cutoff", type=float, default=1, help="Empirical p-value cutoff for including hits")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


if __name__ == '__main__': main()