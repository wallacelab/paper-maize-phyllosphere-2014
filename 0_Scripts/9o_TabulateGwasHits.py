__author__ = 'jgwall'

import argparse
import numpy as np
import pandas as pd

debug = False


def main():
    args = parse_args()
    print("Calculating the number of good hits for traits in",args.infile,"at p-value cutoff of",args.pval_cutoff)
    data=pd.read_table(args.infile)

    # Subset
    data=data.loc[data['empirical_pval'] <= args.pval_cutoff]

    # Count
    traits, counts = np.unique(data['Trait'], return_counts=True)

    # Output
    output = pd.DataFrame({"trait":traits})
    output['count'] = counts
    output.to_csv(args.outfile, sep='\t')

    # Print summary
    summary, summary_counts = np.unique(output['count'], return_counts=True)
    for mysum, mycount in zip(summary, summary_counts):
        print('\t',mycount,"traits have",mysum,"hits")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="File of clustered SNP hits from 4x_ClusterSnpHits.r (with the .best.txt suffix)")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-p", "--pval-cutoff", type=float, default=1, help="Empirical p-value cutoff for including hits")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


if __name__ == '__main__': main()