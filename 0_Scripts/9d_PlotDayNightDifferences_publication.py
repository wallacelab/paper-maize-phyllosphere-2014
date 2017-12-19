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

    # # Write out new table
    # outfile=args.outprefix + ".txt"
    # print("Writing subsetted OTU table with metadata to", outfile)
    # core.to_csv(outfile, sep='\t')

    # Subset to just the targets
    targets = set(args.lines + ['B73'])
    is_target = [g in targets for g in core['genotype']]
    toplot = core.loc[is_target]

    # Remove anything that doesn't have a complete pair (mostly an issue for B73)
    key = np.array([plot + date for plot, date in zip(toplot['plot'], toplot['date'])])
    tally = np.array([np.sum(key == k) for k in key])
    toplot = toplot.loc[tally==2]

    # Plot pairs
    name_key = load_name_key(args.name_key)
    plot_barcharts(toplot, args.outprefix, args.normalize_level, len(args.lines), args.plot_order, name_key)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-c", "--core-otus", help='File with list of core OTUs (Output from step 9c)')
    parser.add_argument("-o", "--outprefix", help='Output file prefix')
    parser.add_argument("-k", "--keyfile", help='QIIME-formatted keyfile')
    parser.add_argument("-n", "--normalize-level", default=10000, help='Number to divide OTU counts by to normalize them (usually the max to get %)')
    parser.add_argument("-p", "--plot-order", default=False, action="store_true", help='Order by field plot instead of by genotype name')
    parser.add_argument("-l", "--lines", nargs="*", help='Which (non-B73) lines to plot. (B73 always by default)')
    parser.add_argument("--name-key", help='File with key to how to display core OTU names in the legend')
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

def plot_barcharts(core, outprefix, normalize, n_lines, plot_order=False, name_key=None):
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
    lines.remove("B73")

    # Make colormap
    colors = cm.Set1(np.linspace(0, 1, n_otus))  # Colormap for the different OTUs
    colors = {o:c for o,c in zip(otus, colors)}

    # Set up figure
    ncol = n_lines
    nrow = 2
    fig = plt.figure(figsize=(3*ncol, 4*nrow))
    grid = gridspec.GridSpec(nrows=nrow, ncols=ncol, hspace=0.4, wspace=0.4)

    # Add barplots for individual plots
    col=0
    bar_width=0.4
    for l in lines:
        ax = fig.add_subplot(grid[0, col], title=l)
        plot_bars(ax, core, otus, genos, l, bar_width, colors)
        col+=1

    # Add B73
    l="B73"
    ax = fig.add_subplot(grid[1, :3], title=l)
    plot_bars(ax, core, otus, genos, l, bar_width, colors)

    # Add legend to colors
    proxies = [patches.Patch(color=colors[o], label=o + "\n(" + name_key[o] + ")") for o in otus.columns]
    # proxies = [patches.Patch(color=colors[o], label="$bf{test}$") for o in otus.columns]
    ax_key = fig.add_subplot(grid[1, 3:])
    legend_properties = {'weight': 'bold', 'size':'xx-small'}#, 'verticalalignment':'center'}
    legend = ax_key.legend(handles=proxies, frameon=False, loc="upper center", prop=legend_properties, ncol=2, mode="expand",
                  borderpad=0, borderaxespad=1, labelspacing=0.9, handleheight=3.25)
    # legend.set_title("OTU Key", prop={'weight':'bold'})
    ax_key.set_title("OTU Key", fontsize="large", fontweight="bold")
    ax_key.axis('off')

    # Save
    fig.savefig(outprefix + ".png", dpi=100)
    fig.savefig(outprefix + ".svg", dpi=100)

def plot_bars(ax, core, otus, genos, myline, bar_width, colors):
    targets = otus.loc[genos == myline, :]
    targetdata = core.loc[genos == myline, :]
    bottom = np.zeros(len(targets))
    xvals = range(len(targets))

    for o in otus.columns:
        myvals = targets[o]
        ax.bar(left=xvals, height=myvals, width=bar_width, bottom=bottom, color=colors[o], align='center')
        bottom += myvals

    # Prettification of axes
    prettify_axes(ax, xvals, targets, targetdata)

def prettify_axes(ax, xvals, targets, targetdata):
    pad = 0.5
    ax.set_xlim([min(xvals) - pad, max(xvals) + pad])
    ax.set_ylim([0, 0.8])
    ax.set_xticks(xvals)
    returns = np.array(['\n'] * len(targets))
    ax.set_xticklabels(targetdata['plot'] + returns + targetdata['time'],
                       fontsize='xx-small', fontweight='bold')
    ax.margins(ymargins=2)
    ax.set_title(ax.get_title(), fontsize="x-large", fontweight="bold")
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('left')

    yticks = [0, 0.2, 0.4, 0.6, 0.8]
    yticklabels=[str(y) for y in yticks]
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize='small', weight='bold')

def load_name_key(infile):
    key = pd.read_table(infile)
    name_key={str(otu):name for otu, name in zip(key['otu'], key['name'])}
    return name_key

if __name__ == '__main__': main()