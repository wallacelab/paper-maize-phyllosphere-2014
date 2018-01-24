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
    clusters=pd.read_table(args.infile)
    # key = pd.read_table(args.key, index_col=0, header=0)

    # Correlate and cluster traits
    for i in range(len(clusters)):

        # Extract data from DataFrame
        mywin = clusters['window'].iloc[i]
        numtraits = clusters['n_traits'].iloc[i]
        score = clusters['score'].iloc[i]
        traits = clusters['traits'].iloc[i].split(',')
        print("\tAnalyzing cluster at", mywin, "with score", score, "from", numtraits,"traits")
        if numtraits != len(traits):
            print("\tWARNING! Number of named traits (", len(traits), ") does not equal number given (", numtraits, ")")

        # Run analyses
        myprefix= args.outprefix + "." + mywin
        plot_correlations(blups, traits, myprefix)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Input file (*.hitcounts.txt from step 5g)")
    parser.add_argument("-o", "--outprefix", help="Output file prefix")
    parser.add_argument("-b", "--blups", help="TASSEL-formatted file of trait BLUPs (to look for correlations)")
    # parser.add_argument("-k", "--key", help="3-column key of metagenome functions (from step 5a)")
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


if __name__ == '__main__': main()