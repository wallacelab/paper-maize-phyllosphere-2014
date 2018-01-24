__author__ = 'jgw87'
"""
Plot heritabilities for publication
"""

import argparse
import gzip
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import re

debug = False

def main():
    args = parse_args()
    print("Plotting metagenome heritability")

    # Load data
    nrows = 100 if debug else None  # Load just a bit if debugging
    heritability = pd.read_table(args.herits, nrows=nrows)
    perms = pd.read_table(args.perms, nrows=nrows).to_dict("list")

    # Filter based on p-value cutoff
    heritability = heritability.loc[heritability['empirical_pval']<=args.pval_cutoff,:]

    # Identify traits with heritability too high (presumably b/c are artifacts)
    too_high = heritability['h2'] > args.max_h2
    if np.sum(too_high) == 0 :
        print("\tWARNING! Supposed to flag artifacts but found", np.sum(too_high), "traits with h2 above",args.max_h2)
    artifacts = heritability.loc[too_high,].copy()
    print("\tIdentified",len(artifacts),"probable artifacts to plot")

    # Grab the top n non-artifact traits
    goodones = heritability.loc[~too_high,].copy()
    goodones = goodones.sort(columns=['empirical_pval', 'h2'], ascending=[True, False])
    goodones = goodones.iloc[:args.num_goodtraits, :]
    print("\tGrabbed the top",len(goodones),"non-artifact traits to plot in comparison")

    # Make data frame of things to plot, flagging if is a likely artifact or not
    artifacts['artifact'] = True
    goodones['artifact'] = False
    toplot = pd.concat([artifacts, goodones])

    # Pick 4 of each to graph
    np.random.seed(args.seed)
    toplot_artifacts = np.sort(np.random.choice(range(len(artifacts)), size=4, replace=False))
    toplot_goodtraits = np.sort(np.random.choice(range(len(goodones)), size=4, replace=False))
    dist_data = pd.concat([artifacts.iloc[toplot_artifacts,:],  goodones.iloc[toplot_goodtraits,:]])
    dist_traits = dist_data['trait']

    # Load blups for histrograms
    blups = pd.read_table(args.blups, skiprows=2, index_col=0)
    blups = blups.loc[:,dist_traits]
    if len(blups.columns) != len(dist_traits):
        print("\tWARNING! Unable to find some traits to plot histograms of:", set(dist_traits) - set(blups.columns))

    # Load key
    tempkey = pd.read_table(args.keyfile)
    key=dict()
    for id, descr in zip(tempkey['id'], tempkey['description']):
        key[id] = descr

    if args.graphfile:
        print("Plotting heritabilities and distributions")
        fig = plt.figure(figsize=args.figsize)  # Scale plot with amount of data

        ax_violins = fig.add_axes([0.18, 0.1, 0.4, 0.8])
        hist_grid = gridspec.GridSpec(nrows=4, ncols=2, left=0.65, hspace=0.5)
        plot_heritabilities(toplot, perms, ax_violins, key)
        plot_distributions(blups, fig, hist_grid, key, dist_data)

        fig.savefig(args.graphfile + ".png", dpi=150)
        fig.savefig(args.graphfile + ".svg", dpi=150)
        plt.close('all')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--herits", help="Pre-processed heritability files from 9k_AssembleMetagenomeHeritData.py")
    parser.add_argument("--perms", help="Pre-processed permutation files from 9k_AssembleMetagenomeHeritData.py")
    parser.add_argument("-b", "--blups", help="BLUPs to be used for plotting distributions of traits")
    parser.add_argument("-g", "--graphfile", help="Output graphic file prefix")
    parser.add_argument("-c", "--pval-cutoff", type=float, default=0, help="Empirical p-value cutoff. Only traits with an empirical p-value of this"
                                                    "or better will be output to the file specified by --passfile")
    parser.add_argument("-f", "--figsize", default=[8, 4], type=float, nargs=2, help="Figure size for output")
    parser.add_argument("-n", "--num-goodtraits", type=int, help="Number of top non-artifact traits to plot")
    parser.add_argument("-k", "--keyfile", help="Metagenome keyfile to give more description of raw COG and KO terms")
    parser.add_argument("-m", "--max-h2", type=float, default=1, help="Filter out traits with h2 above this (b/c thought to be artifacts)")
    parser.add_argument("-s", "--seed", type=int, default=1, help="Random seed for which traits to plot histograms of")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def plot_heritabilities(data, perms, ax, key):
    xleft, xright = -0.05, 1.35
    good_herit = ~data['artifact']
    colors = ["firebrick" if good else "gray" for good in good_herit]

    # Individual scatter plot
    traits, h2 = np.array(data['trait']), np.array(data['h2'])

    ypos = np.arange(len(data))
    markers = ax.scatter(data['h2'], ypos, color=colors, zorder=20, s=60)  # , "o") #, color=np.array(color))

    # Add violinplots for permutations
    violins = list()  # To save violins for later prettification
    for y in ypos:
        if traits[y] in perms:
            violin = ax.violinplot(perms[traits[y]], positions=[y], vert=False)
            violins.append(violin)

    # Add empirical p-values
    text_x = 1.17
    ax.text(x=text_x - 0.01, y=-1, s="p-value", fontsize="xx-small", fontweight="bold", va="center")
    for yval, pval, color in zip(ypos, data['empirical_pval'], colors):
        pval = "{:.4f}".format(pval)
        if color != "gray": color = "black"  # Don't want to colorize same as scatter points
        ax.text(text_x, yval, s=pval, fontsize="xx-small", fontweight="bold", va="center", color=color)

    # Add exact heritability
    text_x = 1.02
    ax.text(x=text_x + 0.04, y=-1, s="h$^2$", fontsize="xx-small", fontweight="bold", va="center")
    for yval, h2, color in zip(ypos, data['h2'], colors):
        if color != "gray": color = "black"  # Don't want to colorize same as scatter points
        h2 = "{:.3f}".format(round(h2, 3))
        ax.text(text_x, yval, s=h2, fontsize="xx-small", fontweight="bold", va="center", color=color)

    # Prettify axes
    ax.set_yticks(ypos)
    ax.set_ylim(bottom=ypos[0] - 1.5, top=ypos[-1] + 1)
    ax.set_xlim(left=xleft, right=xright)
    ax.yaxis.set_ticks_position('none')  # bottom, top, both, or none
    ax.xaxis.set_ticks_position('bottom')  # left, right, both, or none
    ax.invert_yaxis()
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])

    # Prettify X & Y labels
    labels = prettify_labels(traits, key)
    ax.set_yticklabels(labels=labels, fontsize="xx-small", weight="bold")
    ax.set_xlabel(xlabel="Heritability (h$^2$)", fontweight='bold')
    for xlabel in ax.get_xticklabels():
        xlabel.set_fontsize("xx-small")
        xlabel.set_fontweight("bold")

    # Prettify violins
    for v in violins:
        for b in v['bodies']:
            b.set_color("blue")
        for line in [v['cbars'], v["cmins"], v["cmaxes"]]:
            # line.set_color("darkblue")
            line.set_visible(False)

def plot_distributions(blups, fig, grid, key, dist_data):
    nrow, ncol = grid.get_geometry()
    i=0
    for row in range(nrow):
        for col in range(ncol):
            ax = fig.add_subplot(grid[row, col])

            trait = prettify_labels([blups.columns[i]], key=key)[0]
            values = blups.iloc[:,i]
            values = values[np.isfinite(values)]
            isgood = ~dist_data['artifact'].iloc[i]
            n, bins, patches = ax.hist(values, bins=30, color="firebrick" if isgood else "gray")

            # Format axes

            ax.set_title(trait, weight='bold', fontsize='xx-small')
            ax.yaxis.set_ticks_position('none')  # bottom, top, both, or none
            ax.xaxis.set_ticks_position('bottom')  # left, right, both, or none
            ax.tick_params(labeltop='off', labelright='off', labelleft='off', labelsize='xx-small')
            ax.set_ylim(top = ax.get_ylim()[1] * 1.05)  # So bars don't touch the upper line
            ax.xaxis.set_major_locator(plt.MaxNLocator(4))  # Set number of x-ticks

            # Format bars
            for p in patches:
                p.set_linewidth(0.1)

            i+=1


def prettify_labels(labels, key):
    for i in range(len(labels)):

        # flowering time control
        if labels[i] == "flowering_time": labels[i] = "Flowering time"

        # Remove unqiue ID tag and any "log_" prefix
        labels[i] = re.sub(pattern="_[0-9]+$", string=labels[i], repl="")
        labels[i] = re.sub(pattern="^log_", string=labels[i], repl="(log) ")
        labels[i] = re.sub(pattern="_", string=labels[i], repl=" ")

        # # Expand labels for raw KO and COG terms
        # if labels[i] in key:
        #     if key[labels[i]] == "None": key[labels[i]] = "Unknown"
        #     labels[i] = labels[i] + " (" + key[labels[i]] + ")"

        labels[i] = re.sub(pattern="___", string=labels[i], repl = " - ") # fix one formatting issue that was introduced in process of making traits

    return labels


if __name__ == "__main__":
    main()