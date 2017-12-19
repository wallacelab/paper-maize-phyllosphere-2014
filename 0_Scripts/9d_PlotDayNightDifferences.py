__author__ = 'jgwall'

import argparse
import math
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

debug = False
debug_plots=10

def main():
    # Read in data
    args = parse_args()
    otus = pd.read_table(args.infile, index_col=0, low_memory=False)
    core_list = pd.read_table(args.core_otus)
    key = pd.read_table(args.keyfile)

    # Subset to just core OTU table and add metadata from sample name
    core = get_core_otus(otus, core_list)
    core = add_metadata(core, key)

    # Write out new table
    outfile=args.outprefix + ".txt"
    print("Writing subsetted OTU table with metadata to", outfile)
    core.to_csv(outfile, sep='\t')

    # Plot pairs
    plot_barcharts(core, args.outprefix, args.normalize_level, args.plot_order)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-c", "--core-otus", help='File with list of core OTUs (Output from step 9c)')
    parser.add_argument("-o", "--outprefix", help='Output file prefix')
    parser.add_argument("-k", "--keyfile", help='QIIME-formatted keyfile')
    parser.add_argument("-n", "--normalize-level", default=10000, help='Number to divide OTU counts by to normalize them (usually the max to get %)')
    parser.add_argument("-p", "--plot-order", default=False, action="store_true", help='Order by field plot instead of by genotype name')
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def get_core_otus(otus, core_list):
    # Convert OTU names to strings to be clear
    core_names = [str(o) for o in core_list['otu']]
    all_names = [str(o) for o in otus.index]

    # Subset OTU table to just core OTUs
    tokeep = [o in core_names for o in all_names]
    core = otus.loc[tokeep, :]
    core = core.transpose()
    return core

def add_metadata(core, key):
    # Parse QIIMe key to an easy dictionary
    rilkey = {s:str(d) for s,d in zip(key['#SampleID'], key['Description'])}

    # Go through 1 sample at a time and build metadata
    samples = core.index
    time, date, plot, ril = list(), list(), list(), list()
    for s in samples:
        tokens = re.match(string=s, pattern="(LMA.)\.(.+)\.(.......)")
        # Parse sample parts
        mytime = 'unknown'
        if tokens.group(1) == 'LMAD': mytime='day'
        elif tokens.group(1) == 'LMAN': mytime='night'
        mydate = "Aug" + tokens.group(2)
        myplot = tokens.group(3)
        myril = rilkey[s]

        time.append(mytime)
        date.append(mydate)
        plot.append(myplot)
        ril.append(myril)

    # Add to data frame
    core['plot'] = plot
    core['date'] = date
    core['time']=time
    core['genotype'] = ril

    # Sort
    core = core.sort(columns=['genotype','plot','date','time'])
    return core

def plot_barcharts(core, outprefix, normalize, plot_order=False):
    genos = np.array(core['genotype'])
    otus = core.drop(['genotype','plot','date','time'], axis=1) # Just the OTU counts for plotting
    otus/=normalize
    n_otus = len(otus.columns)

    # Determine lines.
    lines=list()
    if plot_order:  # Ordered according to plots in the field
        tmp = core.sort(['plot','date','time'])
        for l in tmp['genotype']:
            if l in lines: continue
            lines.append(l)
    else:   # Ordered alphabetically by line
        lines = sorted(set(genos))

    # Make colormap
    colors = cm.Set1(np.linspace(0, 1, n_otus))  # Colormap for the different OTUs
    colors = {o:c for o,c in zip(otus, colors)}

    # Determine number of plots
    tmp, geno_counts = np.unique(genos, return_counts=True)
    nplots = sum(geno_counts >=2)
    if debug: nplots = debug_plots
    ncol=8
    nrow = math.ceil(nplots / ncol)
    ncol+=1 # Add 1 for legend

    # Set up figure
    fig = plt.figure(figsize=(3*ncol, 4*nrow))
    grid = gridspec.GridSpec(nrows=nrow, ncols=ncol, hspace=0.4, wspace=0.4)

    # Add barplots
    row, col, total = 0,0,0
    bar_width=0.4
    for l in lines:
        count = sum(genos==l)
        if count <2 : continue  # No need to plot a single barchart; no way to compare it
        total +=1
        if debug and total > debug_plots: break
        targets = otus.loc[genos==l,:]
        targetdata = core.loc[genos==l,:]
        ax = fig.add_subplot(grid[row, col], title=l)
        bottom = np.zeros(len(targets))
        xvals = range(len(targets))

        for o in otus.columns:
            myvals = targets[o]
            ax.bar(left = xvals, height=myvals, width=bar_width, bottom=bottom, color=colors[o], align='center')
            bottom+=myvals

        # Prettification of axes
        pad=0.5
        ax.set_xlim([min(xvals)-pad, max(xvals)+pad])
        ax.set_xticks(xvals)
        returns = np.array(['\n'] * len(targets))
        ax.set_xticklabels(targetdata['plot'] + returns + targetdata['date'] + returns + targetdata['time'], fontsize='xx-small', fontweight='bold')
        ax.margins(ymargins=2)


        # Change where next plot will appear (so don't override)
        col+=1
        if col >= ncol-1:
            col = 0
            row +=1

    # Add key to colors
    proxies = [patches.Patch(color=colors[o], label=o) for o in otus.columns]


    ax_key = fig.add_subplot(grid[:, ncol-1], title="Color Key")
    ax_key.legend(handles=proxies, fontsize='x-large')
    ax_key.axis('off')


    # Save
    fig.savefig(outprefix + ".png", dpi=100)
    fig.savefig(outprefix + ".svg", dpi=100)

if __name__ == '__main__': main()