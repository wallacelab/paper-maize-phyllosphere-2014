__author__ = 'jgwall'

import argparse
import math
import numpy as np
import pandas as pd
import re

debug = False
# matplotlib.rcParams.update({"font.size":'10'})

def main():
    args = parse_args()
    print("Getting nearest candidate gene to each hit")

    # Load data
    data = pd.read_table(args.infile)
    flowering_genes = pd.read_csv(args.flowering_genes)

    # Subset hits to just those with right empirical p-value
    hits = data.loc[data['empirical_pval'] <= args.pval_cutoff].copy()  # Copy() to avoid Pandas complaining later
    print("\t",len(hits),"hits have empirical p-values <=",args.pval_cutoff)

    # Subset candidates to make them easier to work with
    candidates = flowering_genes[['gene_name','AGPv3_chrom','AGPv3_start','AGPv3_end']].copy()
    candidates.columns=['gene','chrom','start','stop']
    candidates['pos'] = np.floor((candidates['start'] + candidates['stop'])/2)   # Get an average position

    # Find the nearest candidate for each
    nearest, dist = list(), list()
    for i in range(len(hits)):
        mycandidates = candidates.loc[candidates['chrom'] == hits['Chr'].iloc[i]].copy()
        mydists = np.array(mycandidates['pos'] - hits['Pos'].iloc[i])
        closest = np.argmin(np.abs(mydists))
        nearest.append(mycandidates['gene'].iloc[closest])
        dist.append(mydists[closest])
    hits['nearest_gene'] = nearest
    hits['gene_distance'] = dist

    # Write out
    hits.to_csv(args.outfile, sep='\t', index=None)




def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Input files of collected GWAS results for flowering time")
    parser.add_argument("-o", "--outfile", help="Output file")
    parser.add_argument("-f", "--flowering-genes", help="File of flowering candidate genes")
    parser.add_argument("-p", "--pval-cutoff", default=0.1, type=float, help="Cutoff for empirical p-values to be a 'hit'")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


if __name__ == '__main__': main()