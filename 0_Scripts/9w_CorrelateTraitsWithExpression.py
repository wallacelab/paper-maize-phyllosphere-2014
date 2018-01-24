__author__ = 'jgwall'

import argparse
import pandas as pd
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from seaborn import regplot


debug = False


def main():
    args = parse_args()

    # Load data
    print("Loading data")
    mydata = pd.read_table(args.infile, index_col=0)
    expression = pd.read_table(args.expression)
    key = load_key(args.key)
    print("\tLoaded trait data on total",len(mydata),"traits in",len(mydata.columns),"samples")
    print("\tLoaded expression data on total", len(expression.columns), "genes in", len(expression), "samples")


    # Flip my data so samples are in rows and then filter down
    print("Filtereing to chosen traits and genes")
    mydata = mydata.transpose()
    mydata = mydata.loc[:,args.traits].copy()
    expression = expression.loc[:,args.genes].copy()
    print("\tFiltered to",len(mydata.columns),"traits and", len(expression.columns),"genes")
    if len(mydata.columns) != len(args.traits): print("\tWarning: Some traits were not found")
    if len(expression.columns) != len(args.genes): print("\tWarning: Some genes were not found")

    # Match samples up
    plotdata = match_samples(mydata, expression, key)

    # Set up figure
    print("Plotting correlations")
    nrow = len(args.traits)
    ncol = len(args.genes)
    fig  = plt.figure(figsize=(3*ncol, 3*nrow))
    grid=gridspec.GridSpec(nrows=nrow, ncols=ncol, hspace=0.5, wspace=0.5)
    matplotlib.rcParams.update({'font.size': 8})

    # Plot correlations:
    row, col = 0, 0
    for trait in args.traits:
        for gene in args.genes:
            ax = fig.add_subplot(grid[row, col])

            # build dataset
            xvals, yvals = list(), list()
            for sample in sorted(plotdata.keys()):
                if gene in plotdata[sample] and trait in plotdata[sample]:
                    xvals.append(plotdata[sample][gene])
                    yvals.append(plotdata[sample][trait])
            xvals, yvals = np.array(xvals), np.array(yvals)

            # Get correlation
            isna = np.isnan(xvals) | np.isnan(yvals)
            cor, pval = stats.pearsonr(xvals[~isna], yvals[~isna])
            cor=cor*cor # Get r^2 value

            regplot(x=xvals, y=yvals, scatter=True, ax=ax, color='darkred', scatter_kws={'alpha':0.5, 'linewidths':0, 'color':'blue'})
            ax.text(x=0.05, y=0.95, s="r$^2$ = " + str(cor), fontdict={'fontsize':'x-small', 'weight':'bold'}, transform=ax.transAxes)
            ax.text(x=0.05, y=0.90, s="p = " + str(pval), fontdict={'fontsize': 'x-small', 'weight': 'bold'},
                    transform=ax.transAxes)


            # Prettify axes
            ax.tick_params(labelsize='xx-small', pad=2)
            ax.set_xlabel(gene, fontsize='x-small', weight='bold')
            ax.set_ylabel(trait, fontsize='x-small', weight='bold')

            col+=1
        col=0
        row+=1

    fig.savefig(args.outfile, dpi=150)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Table of trait values from bacterial data")
    parser.add_argument("-e", "--expression", help="Table of RNA expression values")
    parser.add_argument("-k", "--key", help="Keyfile to link expression samples to trait samples")
    parser.add_argument("-t", "--traits", nargs="*", help="Which traits to look at")
    parser.add_argument("-g", "--genes", nargs="*", help="Gene names to look for expression in")
    parser.add_argument("-o", "--outfile", help="Output file (graphical)")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def load_key(infile):
    raw = pd.read_table(infile)
    key=dict()
    for sample, id in zip(raw['expression'], raw['id']):
        if sample in key: print("\tWarning!", sample, "already in key with value", key[sample])
        key[sample]=id
    return(key)


def match_samples(mydata, expression, key):
    print("Matching samples across datasets")

    # I kept trying to get pandas to merge the dataframes but it wouldn't work, possibly due to some samples having multiple gene expression values somehow
    expression['id'] = [key[s] for s in expression.index]
    sampledata=dict()

    # Add trait data
    for sample in mydata.index:
        if sample not in sampledata: sampledata[sample]=dict()
        for trait in mydata.columns:
            if trait in sampledata[sample]:
                print("\tWARNING! Trait",trait,"already loaded for sample",sample)
            sampledata[sample][trait] = mydata.loc[sample, trait]

    # Add gene expression data
    for i in range(len(expression)):
        sample = expression['id'].iloc[i]
        if sample not in sampledata: sampledata[sample] = dict()

        for gene in expression.columns:
            if gene == 'id' : continue
            if gene in sampledata[sample]:
                print("\tWARNING! Gene",gene,"already loaded for sample",sample," and will be set to NaN")
                sampledata[sample][gene] = np.nan
            else:
                sampledata[sample][gene] = expression[gene].iloc[i]

    return sampledata


# def prettify_axes(ax):
#
#     # Shrink and prettify tick labels
#     xticklabels = list()
#     for tick in ax.get_xticks():
#         xticklabels.append(str(tick))
#     ax.set_xticklabels(xticklabels, weight='bold', size='small')
#     # y-axis
#     yticklabels = list()
#     for tick in ax.get_yticks():
#         yticklabels.append(str(tick))
#     ax.set_yticklabels(yticklabels, weight='bold', size='small')

if __name__ == '__main__': main()