__author__ = 'jgwall'

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
from sys import exit
from biokit.viz import corrplot

debug = False


def main():
    args = parse_args()
    print("Analyzing clusters of hits in",args.infile)

    # Load data
    blups=pd.read_table(args.blups, index_col=0, skiprows=2)
    scores=pd.read_table(args.infile)
    key = pd.read_table(args.key, index_col=0, header=0)

    # Filter to just target clusters
    tokeep = (scores['cutoff']==args.pval_cutoff) & (scores['corrected_count'] >= args.min_score)
    clusters = scores.loc[tokeep]
    print("After filtering for empirical p=",args.pval_cutoff,"and corrected count >=",args.min_score,"have",len(clusters),"clusters remaining")
    if len(clusters) == 0:
        print("\tERROR: NO CLUSTERS FOUND. Cannot perform analysis so script is exiting")
        exit()

    # Correlate and cluster traits
    outprefix = args.outprefix + ".p" + str(args.pval_cutoff) + ".score" + str(args.min_score)
    for i in range(len(clusters)):

        # Extract data from DataFrame
        mychrom = clusters['chrom'].iloc[i]
        mypos = int(clusters['window_pos'].iloc[i])
        numtraits = clusters['num_traits'].iloc[i]
        score = clusters['corrected_count'].iloc[i]
        traits = clusters['traits'].iloc[i].split(',')
        print("\tAnalyzing cluster at", str(mychrom) + ":" + str(mypos), "with score", score, "from", numtraits,
              "traits")
        if numtraits != len(traits):
            print("\tWARNING! Number of named traits (", len(traits), ") does not equal number given (", numtraits, ")")

        # Run analyses
        myprefix= outprefix + "." + str(mychrom) + "_" + str(mypos)
        plot_correlations(blups, traits, myprefix)
        parse_annotations(traits, key,  myprefix)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Input file (*.hitcounts.txt from step 5g)")
    parser.add_argument("-o", "--outprefix", help="Output file prefix")
    parser.add_argument("-b", "--blups", help="TASSEL-formatted file of trait BLUPs (to look for correlations)")
    parser.add_argument("-k", "--key", help="3-column key of metagenome functions (from step 5a)")
    parser.add_argument("-s", "--min-score", type=float, default=1, help="Minimum score for a SNP cluster to be analyzed")
    parser.add_argument("-p", "--pval-cutoff", type=float, default=0.05, help="Empirical p-value cutoff to use")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def plot_correlations(blups, traits, outprefix):
    # Subset blups to just the traits listed
    myblups = blups[traits]

    # Set up figure
    outpng = outprefix + ".corrplot.png"
    fig = plt.figure(figsize=(len(traits)*1.25, len(traits)))
    ax=fig.add_subplot(111)

    # Plot correlation matrix
    cors = corrplot.Corrplot(myblups)
    cors.plot(ax=ax, lower='ellipse', upper='number')
    fig.savefig(outpng, dpi=100)


def parse_annotations(traits, key, outprefix):
    mytraits = [re.sub(string=t, pattern="^log_", repl="") for t in traits]
    mytraits = [re.sub(string=t, pattern="_[0-9]+$", repl="") for t in mytraits]
    mykey = key.loc[mytraits,:]
    mykey.to_csv(outprefix + ".annotations.txt", sep='\t', index_label="id")


if __name__ == '__main__': main()