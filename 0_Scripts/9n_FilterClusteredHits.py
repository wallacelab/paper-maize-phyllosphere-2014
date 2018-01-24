__author__ = 'jgwall'

import argparse
import numpy as np
import pandas as pd

debug = False


def main():
    args = parse_args()
    print("Filtereing GWAS hits in",args.infile,"based on heritability (to exclude artifacts) in",args.herits)

    # Load data
    gwas = pd.read_table(args.infile)
    herits = pd.read_table(args.herits)

    # Get list of bad traits
    badtraits = set(herits['trait'].loc[herits['h2']>= args.max_herit])

    # Filter gwas clusters
    badhit = np.array([mytrait in badtraits for mytrait in gwas['Trait']])
    filtered = gwas.loc[~badhit, :]
    removed = set(gwas['Trait'].loc[badhit])
    print("\tRemoved",sum(badhit),"hits across",len(removed),"traits (out of",len(gwas),"total hits)")

    # Write out
    filtered.to_csv(args.outfile, sep='\t', index=False)
    badtraits = sorted(badtraits)
    OUT = open(args.traitfile, "w")
    [OUT.write(t + "\n") for t in badtraits]
    OUT.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Clustered GWAS hits")
    parser.add_argument("--herits", help="Heritability file from step 9k")
    parser.add_argument("-o", "--outfile", help="Output file of filetered cluster results")
    parser.add_argument("-t", "--traitfile", help="Output file of 'bad' traits that were excluded based on filtering")
    # parser.add_argument("-p", "--max-pval", type=float, default=1, help="Maximum p-value to be considered a good hit")
    parser.add_argument("-m", "--max-herit", type=float, default=1, help="Exclude heritability above this (b/c is probably a artifact)")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()




if __name__ == '__main__': main()