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

    # Filter out traits with heritability too high (presumably b/c are artifacts)
    too_high = heritability['h2'] > args.max_h2
    if np.sum(too_high) > 0 :
        print("\tFiltering out", np.sum(too_high), "traits with h2 above",args.max_h2,"because assumed to be artifacts" )
        heritability = heritability.loc[~too_high,]

    # Filter out non-target traits
    if args.traits:
        heritability = subset_to_traits(heritability, args.traits)

    # Load key
    tempkey = pd.read_table(args.keyfile)
    key=dict()
    for id, descr in zip(tempkey['id'], tempkey['description']):
        key[id] = descr

    # Determine how many traits pass p-value filter
    good_herit = heritability['empirical_pval'] <= args.pval_cutoff
    print("\t",np.sum(good_herit),"traits have empirical p-value of",args.pval_cutoff,"or better")

    if args.num_traits:
        print("Reducing to just top",args.num_traits,"traits")
        heritability = heritability.iloc[:args.num_traits, :]

    if args.graphfile:
        plot_heritabilities(heritability, perms, args.graphfile, args.pval_cutoff, args.figsize, key)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--herits", help="Pre-processed heritability files from 9k_AssembleMetagenomeHeritData.py")
    parser.add_argument("--perms", help="Pre-processed permutation files from 9k_AssembleMetagenomeHeritData.py")
    parser.add_argument("-g", "--graphfile", help="Output graphic file prefix")
    parser.add_argument("-c", "--pval-cutoff", type=float, default=0, help="Empirical p-value cutoff. Only traits with an empirical p-value of this"
                                                    "or better will be output to the file specified by --passfile")
    parser.add_argument("-f", "--figsize", default=[8, 4], type=float, nargs=2, help="Figure size for output")
    parser.add_argument("-n", "--num-traits", type=int, help="Number of top traits to plot")
    parser.add_argument("-k", "--keyfile", help="Metagenome keyfile to give more description of raw COG and KO terms")
    parser.add_argument("-m", "--max-h2", type=float, default=1, help="Filter out traits with h2 above this (b/c thought to be artifacts)")
    parser.add_argument("-t", "--traits", help="File of traits to filter for")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def plot_heritabilities(data, perms, graphfile, cutoff, figsize, key):
    fig = plt.figure(figsize=figsize)  # Scale plot with amount of data
    xleft, xright = -0.05, 1.2
    good_herit = data['empirical_pval'] <= cutoff
    colors = ["firebrick" if good else "gray" for good in good_herit]

    # Individual scatter plot
    traits, h2 = np.array(data['trait']),np.array(data['h2'])
    ax=fig.add_axes([0.5, 0.1, 0.45, 0.8])
    ypos=np.arange(len(data))
    markers = ax.scatter(data['h2'], ypos, color=colors, zorder=20, s=60)#, "o") #, color=np.array(color))

    # Add violinplots for permutations
    violins = list()    # To save violins for later prettification
    for y in ypos:
        if traits[y] in perms:
            violin = ax.violinplot(perms[traits[y]], positions=[y], vert=False)
            violins.append(violin)

    # Add empirical p-values
    text_x = 1.02
    ax.text(x=text_x, y = -1, s = "p-value", fontsize="xx-small", fontweight="bold", va="center")
    for yval, pval, color in zip(ypos, data['empirical_pval'], colors):
        if color != "gray": color = "black" # Don't want to colorize same as scatter points
        ax.text(text_x, yval, s="{:.4f}".format(pval), fontsize="xx-small", fontweight="bold", va="center", color=color)



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
            #line.set_color("darkblue")
            line.set_visible(False)

    fig.savefig(graphfile + ".png", dpi=150)
    fig.savefig(graphfile + ".svg", dpi=150)
    plt.close('all')

def prettify_labels(labels, key):
    for i in range(len(labels)):

        # flowering time control
        if labels[i] == "flowering_time": labels[i] = "Flowering time"

        # Remove unqiue ID tag and any "log_" prefix
        labels[i] = re.sub(pattern="_[0-9]+$", string=labels[i], repl="")
        labels[i] = re.sub(pattern="^log_", string=labels[i], repl="")

        # Expand labels for raw KO and COG terms
        if labels[i] in key:
            if key[labels[i]] == "None": key[labels[i]] = "Unknown"
            labels[i] = labels[i] + " (" + key[labels[i]] + ")"

        labels[i] = re.sub(pattern="___", string=labels[i], repl = " - ") # fix one formatting issue that was introduced in process of making traits

    return labels

def subset_to_traits(heritability, traitfile):
    print("\tFiltering to the traits in", traitfile)

    # Make set of traits to filter with
    traits = open(traitfile).readlines()
    traits = [t.strip().replace(" ", "_").replace("-", "_") for t in traits]  # Fix formatting
    traits = set(traits)

    # Adjust h2 trait labels so match
    mytraits = [re.sub(string=t, pattern="_[0-9]+$", repl="") for t in heritability['trait']]  # Remove unique tag on end
    mytraits = [re.sub(string=t, pattern="^log_", repl="") for t in mytraits]  # Remove unique tag on end
    tokeep = [t in traits for t in mytraits]

    # Filter
    heritability = heritability.loc[tokeep, :]

    return(heritability)

if __name__ == "__main__":
    main()