__author__ = 'jgw87'
"""
Take a list of TASSEL files with permutations and calculate heritabilities
"""

import argparse
import gzip
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import re


def main():
    args = parse_args()
    heritability, perms = calc_heritability(args.infiles)
    if args.outfile:
        heritability.to_csv(args.outfile, sep='\t', index=False)
    if args.graphfile:
        plot_heritabilities(heritability, perms, args.graphfile, args.pval_cutoff)
    if args.passfile:
        output_good_traits(heritability, args.pval_cutoff, args.passfile)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", help="TASSEL MLM output files", nargs="*")
    parser.add_argument("-o", "--outfile", help="Output text file")
    parser.add_argument("-g", "--graphfile", help="Output graphic")
    parser.add_argument("-c", "--pval-cutoff", type=float, default=0, help="Empirical p-value cutoff. Only traits with an empirical p-value of this"
                                                    "or better will be output to the file specified by --passfile")
    parser.add_argument("-p", "--passfile", help="Output file for the names of traits that pass the filter in --pval-cutoff")
    args = parser.parse_args()
    return args

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

def plot_heritabilities(data, perms, graphfile, cutoff):
    fig = plt.figure(figsize=(4 + 0.15 * len(data),6))  # Scale plot with amount of data
    grid = gridspec.GridSpec(nrows=100, ncols=100)
    vert = 80
    ybottom, ytop = -0.05, 1.05
    colors = ["blue" if pval <= cutoff else "gray" for pval in data['empirical_pval']]

    # Violin plot of distribution
    ax1=fig.add_subplot(grid[:vert,:25], title="Distribution of heritabilities")
    ax1.violinplot(data['h2'])
    ax1.set_xticklabels(labels="", visible=False)
    ax1.set_ylim(bottom=ybottom, top=ytop)

    # Individual plot
    traits, h2 = np.array(data['trait']),np.array(data['h2'])
    ax2=fig.add_subplot(grid[:vert,30:], title="Exact heritabilities", xlabel="trait", ylabel="Heritability")
    xpos=np.arange(len(data))
    markers = ax2.scatter(xpos, data['h2'], color=colors, zorder=20)#, "o") #, color=np.array(color))
    ax2.set_xticks(xpos)
    ax2.set_xticklabels(labels=traits, rotation=90, fontsize="xx-small")
    ax2.set_xlim(left=xpos[0]-1, right=xpos[-1]+1)
    ax2.set_ylim(bottom=ybottom, top=ytop)

    # Add boxplots for permutations
    for x in xpos:
        if traits[x] in perms:
            ax2.violinplot(perms[traits[x]], positions=[x])

    fig.savefig(graphfile, dpi=150)
    plt.close('all')

def output_good_traits(data, cutoff, outfile):
    print("Outputting good traits to",outfile)
    tokeep = (data['empirical_pval'] <= cutoff) & (~pd.isnull(data['h2']))    # Had to add null filter b/c broken analyses still passed
    output = data.loc[tokeep,:]
    output.to_csv(outfile, sep='\t', index=False)


if __name__ == "__main__":
    main()