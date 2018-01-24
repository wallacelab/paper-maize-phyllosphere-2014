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
    if debug: args.infiles = args.infiles[:5]
    heritability, perms = calc_heritability(args.infiles)

    if args.num_traits:
        print("Reducing to just top",args.num_traits,"traits")
        heritability = heritability.iloc[:args.num_traits, :]

    if args.outfile:
        heritability.to_csv(args.outfile, sep='\t', index=False)
    if args.graphfile:
        plot_heritabilities(heritability, perms, args.graphfile, args.pval_cutoff, args.figsize)
    if args.passfile:
        output_good_traits(heritability, args.pval_cutoff, args.passfile)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", help="TASSEL MLM output files", nargs="*")
    parser.add_argument("-o", "--outfile", help="Output text file")
    parser.add_argument("-g", "--graphfile", help="Output graphic file prefix")
    parser.add_argument("-c", "--pval-cutoff", type=float, default=0, help="Empirical p-value cutoff. Only traits with an empirical p-value of this"
                                                    "or better will be output to the file specified by --passfile")
    parser.add_argument("-p", "--passfile", help="Output file for the names of traits that pass the filter in --pval-cutoff")
    parser.add_argument("-f", "--figsize", default=[8, 4], type=float, nargs=2, help="Figure size for output")
    parser.add_argument("-n", "--num-traits", type=int, help="Number of top traits to plot")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def calc_heritability(infiles):
    print("Calculating heritabilities in",len(infiles),"files")

    # First make a master dict() of all heritabilities; will split into perms and originals later
    heritability = dict()
    for infile in infiles:
        #print(infile)
        if len(open(infile).readline()) ==0 : continue #Skip empty files
        stats=pd.read_csv(infile, sep='\t')
        for trait, gen_var, resid_var in zip(stats['Trait'], stats['Genetic Var'], stats['Residual Var']):
            #if 'med1000_otu10_Rarefy10000_unifrac_bdiv_even10000_unweighted_unifrac_PC1' in trait:
                #print(trait,"has Vg", gen_var,"and Vr",resid_var)
            #print(trait)
            if trait in heritability: continue  # Already loaded, so skip
            heritability[trait] = gen_var / (gen_var + resid_var)

    # Now split into original traits and permutations
    core, perms = dict(), dict()
    for trait in sorted(heritability.keys()):
        if "_perm" in trait:
            basename = re.sub(string=trait, pattern="_perm.+", repl="")
            if basename not in perms:
                perms[basename] = list()
            perms[basename].append(heritability[trait])
        else:
            core[trait] = heritability[trait]


    traits = sorted(core.keys())
    for t in traits:
        if t not in perms: perms[t] = [np.nan, np.nan]  # Hack to get around traits with no permutation data
    h2 = [core[t] for t in traits]
    means = [np.nanmean(perms[t]) for t in traits]
    medians = [np.nanmedian(perms[t]) for t in traits]
    maxes = [np.nanmax(perms[t]) for t in traits]
    pvals = [np.sum(perms[t] >= core[t])/len(perms[t]) for t in traits] # Empirical p-value is how often got a heritability that high or more
    heritability = pd.DataFrame({"trait":traits, "h2":h2, "perm_mean":means, "perm_median":medians, "perm_max":maxes, "empirical_pval":pvals})
    heritability = heritability[["trait", "h2", "perm_mean", "perm_median","perm_max","empirical_pval"]]
    return heritability.sort(columns=["empirical_pval", "h2"], ascending=[True, False]), perms

def plot_heritabilities(data, perms, graphfile, cutoff, figsize):
    fig = plt.figure(figsize=figsize)  # Scale plot with amount of data
    xleft, xright = -0.05, 1.35
    colors = ["firebrick" if pval <= cutoff else "gray" for pval in data['empirical_pval']]

    # Individual scatter plot
    traits, h2 = np.array(data['trait']),np.array(data['h2'])
    ax=fig.add_axes([0.3, 0.1, 0.65, 0.8])
    ypos=np.arange(len(data))
    markers = ax.scatter(data['h2'], ypos, color=colors, zorder=20, s=60)#, "o") #, color=np.array(color))

    # Add violinplots for permutations
    violins = list()    # To save violins for later prettification
    for y in ypos:
        if traits[y] in perms:
            violin = ax.violinplot(perms[traits[y]], positions=[y], vert=False)
            violins.append(violin)

    # Add empirical p-values
    text_x = 1.17
    ax.text(x=text_x-0.01, y = -1, s = "p-value", fontsize="xx-small", fontweight="bold", va="center")
    for yval, pval, color in zip(ypos, data['empirical_pval'], colors):
        pval = "{:.4f}".format(pval)
        if color != "gray": color = "black" # Don't want to colorize same as scatter points
        ax.text(text_x, yval, s=pval, fontsize="xx-small", fontweight="bold", va="center", color=color)

    # Add exact heritability
    text_x = 1.02
    ax.text(x=text_x+0.04, y=-1, s="h$^2$", fontsize="xx-small", fontweight="bold", va="center")
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
    labels = prettify_labels(traits)
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

def prettify_labels(labels):
    for i in range(len(labels)):
        # flowering time control
        if labels[i] == "flowering_time": labels[i] = "Flowering time"
        # Families and orders
        elif labels[i].endswith("ales") or labels[i].endswith("aceae"):
            labels[i] = re.sub(pattern=".+\\.", string=labels[i], repl="")
        # OTUs down to species
        labels[i] = re.sub(pattern="(.+)\\.([0-9]+)", string=labels[i], repl="OTU \\2 (\\1)")
    return labels

# # Make my traits look good (italics, etc.)
# def customize_trait(label):
#     text = label.get_text()
#     if text == "flowering_time":
#         label.set_text("Flowering time")
#     elif text.endswith("ales") or text.endswith("aceae"):
#         text = re.sub(pattern = ".+\\.", string=text, repl="")
#         label.set_text(text)


# def plot_heritabilities(data, perms, graphfile, cutoff, figsize):
#     fig = plt.figure(figsize=figsize)  # Scale plot with amount of data
#     vert = 80
#     ybottom, ytop = -0.05, 1.05
#     colors = ["firebrick" if pval <= cutoff else "gray" for pval in data['empirical_pval']]
#
#     # Individual scatter plot
#     traits, h2 = np.array(data['trait']),np.array(data['h2'])
#     ax=fig.add_axes([0.1, 0.3, 0.8, 0.6])
#     xpos=np.arange(len(data))
#     markers = ax.scatter(xpos, data['h2'], color=colors, zorder=20, s=60)#, "o") #, color=np.array(color))
#
#     # Add violinplots for permutations
#     violins = list()    # To save violins for later prettification
#     for x in xpos:
#         if traits[x] in perms:
#             violin = ax.violinplot(perms[traits[x]], positions=[x])
#             violins.append(violin)
#
#
#     # Prettify axes
#     ax.set_xticks(xpos)
#     ax.set_xlim(left=xpos[0] - 1, right=xpos[-1] + 1)
#     ax.set_ylim(bottom=ybottom, top=ytop)
#     ax.xaxis.set_ticks_position('none')  # bottom, top, both, or none
#     ax.yaxis.set_ticks_position('left')  # left, right, both, or none
#
#     # Prettify X & Y labels
#     ax.set_xticklabels(labels=traits, rotation=90, fontsize="xx-small", weight="bold")
#     ax.set_ylabel(ylabel="Heritability (h$^2$)", fontweight='bold')
#     for ylabel in ax.get_yticklabels():
#         ylabel.set_fontsize("xx-small")
#         ylabel.set_fontweight("bold")
#
#     # Prettify violins
#     for v in violins:
#         for b in v['bodies']:
#             b.set_color("blue")
#         for line in [v['cbars'], v["cmins"], v["cmaxes"]]:
#             #line.set_color("darkblue")
#             line.set_visible(False)
#
#     fig.savefig(graphfile, dpi=150)
#     plt.close('all')

def output_good_traits(data, cutoff, outfile):
    print("Outputting good traits to",outfile)
    tokeep = (data['empirical_pval'] <= cutoff) & (~pd.isnull(data['h2']))    # Had to add null filter b/c broken analyses still passed
    output = data.loc[tokeep,:]
    output.to_csv(outfile, sep='\t', index=False)


if __name__ == "__main__":
    main()