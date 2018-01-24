__author__ = 'jgwall'

import argparse
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import re

debug = False
# matplotlib.rcParams.update({"font.size":'10'})

def main():
    args = parse_args()
    data = [pd.read_table(i) for i in args.infiles]
    flowering_genes = pd.read_csv(args.flowering_genes)

    # Set up figure
    fig = plt.figure(figsize=(10, 6.5))
    grid = gridspec.GridSpec(nrows=3, ncols=1, hspace=0.5, wspace=0.2)
    chromlengths = load_chromlengths(args.chromlengths)

    # Plot Manhatten plots
    for i in range(len(data)):
        trait = data[i]
        islast =  i==len(data)-1  # Variable to tell if is last plot for adding legend
        ax = fig.add_subplot(grid[i,:], ylabel="-log$_{10}$ p-value")
        plot_manhatten(ax, trait, chromlengths, add_legend=islast, sparsify=args.sparsify, candidates=flowering_genes, winsize=args.candidate_window)


    fig.savefig(args.outprefix + ".png", dpi=150)
    fig.savefig(args.outprefix + ".svg", dpi=300)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*", help="Input files of collected GWAS results")
    parser.add_argument("-o", "--outprefix", help="Output graphic file prefix")
    parser.add_argument("-c", "--chromlengths", help="File of chromosome lengths")
    parser.add_argument("-f", "--flowering-genes", help="File of flowering candidate genes")
    parser.add_argument("-s", "--sparsify", type=float, default=None, help="What fraction of background GWAS hits to plot")
    parser.add_argument("-w", "--candidate-window", type=int, default=1000000, help="Distance from a hit for a candidate gene to be plotted")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def load_chromlengths(infile):
    lengths=dict()
    for line in open(infile, "r"):
        chrom, length, cum = line.strip().split('\t')[:3]   # ':2' in case there are additional columns for some reas
        if chrom.isnumeric(): chrom = int(chrom)
        lengths[chrom] = {"length":int(length), "cum":int(cum)}
    return lengths


def plot_manhatten(ax, results, chromlengths, sparsify=None, add_legend=False, candidates=None, winsize=1000000):

    # Reduce least significant hits, always keeping the top 100 and ones with empirical p-value flags
    if sparsify:
        newsize = math.floor(sparsify * len(results))
        np.random.seed(1)
        tokeep = np.random.choice(range(len(results)), size=newsize, replace=False)
        top_100 = np.argsort(np.array(results['p']))[:100]

        # Always keep ones that passed empirical p-value threshold, too
        if 'empirical_pval' in results.columns:
            pass_empirical = np.where(results['empirical_pval'] <= 0.1)[0]
        else: pass_empirical = []

        tokeep = sorted(set(tokeep) | set(top_100) | set(pass_empirical))
        results = results.iloc[tokeep,:].copy()

    # Remove NAN values
    results = results.loc[~np.isnan(results['p']),].copy()

    # Colors and size
    colors = 'black'
    size=5
    if 'empirical_pval' in results.columns:
        empiricals = np.array(results['empirical_pval'])   # To speed calculations up
        colors = np.array(['darkgray'] * len(results))
        colors[np.array(results['Chr']) % 2 ==0] = 'gray'    # To delimit chromosomes
        colors[empiricals <= 0.1] = 'black'
        colors[empiricals <= 0.05] = 'blue'
        colors[empiricals <= 0.01] = 'darkred'

        # Adjust size
        size = np.array([5] * len(results))
        size[empiricals <= 0.1] = 20
        size[empiricals <= 0.05] = 30
        size[empiricals <= 0.01] = 40
    results['size'] = size
    results['color'] = colors

    # x-values
    xvals = list()
    for mychrom, mypos in zip(results['Chr'], results['Pos']):
        xvals.append(mypos + chromlengths[mychrom]['cum'])
    results['x'] = xvals

    # Plot background
    background = size == 5
    myscatter = ax.scatter(x=results['x'][background], y=-np.log10(results['p'])[background], s=results['size'][background],
                           color=results['color'][background], alpha=0.25, linewidths=0)
    myscatter.set_rasterized(True)
    # Plot significant hits
    ax.scatter(x=results['x'][~background], y=-np.log10(results['p'])[~background], s=results['size'][~background],
                           color=results['color'][~background], alpha=0.75, linewidths=0)

    # Add chromosome divisions
    for c in chromlengths:
        ax.axvline(chromlengths[c]['cum'], color="lightgray", zorder=-1)
    ax.axvline(chromlengths[c]['cum'] + chromlengths[c]['length'] , color="lightgray", zorder=-1)

    # Adjust x and y limits
    min_x, max_x = min(results['x']), max(results['x'])
    pad = (max_x - min_x)/20
    ax.set_xlim(left=min_x-pad, right=max_x+pad)
    ax.set_ylim(bottom=0)

    # Adjust title
    mytrait = " ".join(np.unique(results['Trait'])) # Join is just so it's obvious if >1 trait got put in somehow
    mytrait = re.sub(string = mytrait, pattern="^log_", repl="log ")
    mytrait = re.sub(string=mytrait, pattern="_[0-9]+$", repl="")   # Strip trailing numbers added during analysis for uniqueness
    # Manual fixes
    if mytrait=="flowering_time": mytrait="Flowering Time"
    mytrait = re.sub(string=mytrait, pattern="K02769", repl="K02769 (PTS system, fructose-specific IIB component)")
    mytrait = re.sub(string=mytrait, pattern="COG0181", repl="COG0181 (Porphobilinogen deaminase)")
    #Set title
    ax.set_title(mytrait, weight="bold", size='medium')

    # Alter axis lines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Alter tick marks
    ax.tick_params(axis='y', which='both', right='off', labelsize='x-small')
    ax.tick_params(axis='x', which='both', top='off', bottom='off')

    # Plot chromosome names
    cnames, cpos, chromdiv = list(), list(), set()
    for mychrom in sorted(chromlengths.keys()):
        if mychrom in {'UNKNOWN', 'Pt', 'Mt', '0'}: continue  # skip inconsequential ones
        cnames.append("chr" + str(mychrom))
        cpos.append(chromlengths[mychrom]['cum'] + chromlengths[mychrom]['length'] / 2)
    ax.set_xticks(cpos)
    ax.set_xticklabels(cnames, fontsize="small")

    # Add legend
    if add_legend:
        dotdummy = patches.Patch(color='green', label='Empirical p-values:', alpha=0)   # Dummy one to put title to left
        dot01 = patches.Patch(color='darkred', label='p ≤ 0.01')
        dot05 = patches.Patch(color='blue', label='p ≤ 0.05')
        dot10 = patches.Patch(color='black', label='p ≤ 0.10')
        ax.legend(handles=[dotdummy, dot01, dot05, dot10], loc='center left', frameon=False,
                           bbox_to_anchor=[-0.05, -0.3], ncol=4, prop={'weight':'bold', 'size':'x-small'})

    # Add candidate genes
    if mytrait == "Flowering Time":
        besthits = results.loc[empiricals <= 0.01, :]
        add_candidate_genes(ax, candidates, besthits, chromlengths, winsize)


def add_candidate_genes(ax, genes, besthits, chromlengths, winsize):
    print("\tAdding candidate genes")
    # subset to be easier to work with
    candidates = genes[['gene_name','AGPv3_chrom','AGPv3_start','AGPv3_end']].copy()
    candidates.columns=['gene','chrom','start','stop']

    # Get plotting position
    candidates['pos'] = (candidates['start'] + candidates['stop'])/2    # Plot at average position
    offsets =  [chromlengths[c]['cum'] for c in candidates['chrom']]
    candidates['xval'] = np.array(offsets) + candidates['pos']

    # Reduce to just candidates that are close to the best hits from plotting
    toplot = [False] * len(candidates)
    target_chroms, target_pos = np.array(besthits['Chr']), np.array(besthits['Pos'])
    for i in range(len(candidates)):
        mychrom = candidates['chrom'].iloc[i]
        mypos = candidates['pos'].iloc[i]
        isnear = (target_chroms == mychrom) & (abs(target_pos - mypos) < winsize)    # Go for candidates within 10 kb on same chromosome
        toplot[i] = np.any(isnear)
    candidates['toplot'] = toplot
    print("\t\tFound",sum(candidates['toplot']),"out of",len(candidates),"candidate genes to plot")
    candidates = candidates.loc[candidates['toplot']]

    # Plot
    for mygene, myx in zip(candidates['gene'], candidates['xval']):
        ax.axvline(myx, color="darkgreen", linewidth=1, zorder=0, alpha=0.75, linestyle='dashed')
        ymax = ax.get_ylim()[1]
        ax.text(x=myx, y=ymax * 0.9, s=mygene, horizontalalignment='left', fontsize='xx-small', weight='bold', color='darkgreen')



def set_color(p, cutoff, chrom):
    if p < cutoff and chrom %2 == 0: return "blue"
    if p < cutoff and chrom %2 != 0: return "red"
    if chrom % 2 == 0: return "gray"
    if chrom % 2 != 0: return "silver"
    return "purple" # Should never get here

if __name__ == '__main__': main()