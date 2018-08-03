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
matplotlib.rcParams.update({"font.size":'10'})

def main():
    args = parse_args()
    results = [pd.read_table(i) for i in args.infiles]
    print("\tLoaded GWAS results from", len(args.infiles), "input files")

    # Remove rows not associated with markers
    for i in range(len(results)):
        r = results[i]
        results[i] = r.loc[ ~np.isnan(r['Chr']),:].copy()   # Copy so don't have issues with setting stuff on slices
        # print("\t\tBefore removing NAN had",len(r),"rows of data; afterward have",len(results[i]))
        if (len(results[i]) - len(r)) > 1:
            print("\t\tWARNING! When removing <NA> chromosones, had", len(results[i]) - len(r), "rows removed (expect 0 or 1)")

    # Confirm only have a single trait
    traits = [np.array(r['Trait']) for r in results]
    traits=np.unique(np.concatenate(traits))  # collapse
    if len(traits) != 1:
        print("\tWARNING!!! Expected 1 trait but got",len(traits),":",traits,"; results may not be interpretable")
    trait="+".join(traits)
    # print(trait)

    # Load simulated QTL
    qtl = load_qtl(args.qtl, args.hapmap, trait)

    # Get adjusted p-values if permutations supplied
    print("\tCalculating empirical p-values based on permutation p-values in",args.perms)
    results = calc_empirical_pvals(results, args.perms)

    # # Filter to just empirical p-values I care about
    outtext = args.outprefix + ".found_qtl.txt"
    print("\t\tWriting out QTL with empirical p-values less than", args.empirical_p_cutoff,"to",outtext)
    output = pd.concat(results)

    # Filter to just QTL
    qtl_names = set(qtl['rs#'])
    is_qtn = [m in qtl_names for m in output['Marker']]
    output = output.loc[is_qtn,:]
    print("\t\t\tFound",len(output),'QTN')

    # Filter to empirical p-values
    output = output.loc[output['empirical_pval'] <= args.empirical_p_cutoff,:]
    if hasattr(pd, "sort"): # Have to test for sort() because got deprecated after pandas 0.17
        output = output.sort(columns=['empirical_pval', 'Chr', 'Pos'])   # Pandas 0.16 and prior
    else:
        output = output.sort_values(by=['empirical_pval', 'Chr', 'Pos'], axis='index')  # Pandas 0.17 and later
    print("\t\t\tFound", len(output), 'QTN with empirical p-values <=',args.empirical_p_cutoff)

    # Print out
    output.to_csv(outtext, sep='\t')

    # Make plot
    results = pd.concat(results)
    chrom_offsets = load_chromlengths(args.chromlengths)
    plot_gwas(results, chrom_offsets, args.outprefix, args.rasterize, args.sparsify, trait, qtl, args.empirical_p_cutoff)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*", help="List of files with GWAS results; requires columns Trait, Chr, Pos, and p")
    parser.add_argument("-o", "--outprefix", help="Output file prefix")
    parser.add_argument("-c", "--chromlengths", help="File of chromosome lengths")
    parser.add_argument("-q", "--qtl", help="File of simulated QTL")
    parser.add_argument("-m", "--hapmap", help="Hapmap file that includes QTL (to get physical positions)")
    parser.add_argument("-p", "--perms", help="File of permuted GWAS results for calculating empirical p-values")
    parser.add_argument("--empirical-p-cutoff", type=float, default=1, help="Where to cut off p-values to output based on their empirical p-value")
    parser.add_argument("-r", "--rasterize", default=True, action="store_true",
                        help="Whether to rasterize the inner portion of each plot. (Only applicable to vector output formats)")
    parser.add_argument("--sparsify", type=float, default=0.1, help="Randomly subset low-significance p-values to this fraction of original to make plotting faster")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def load_chromlengths(infile):
    print("Loading chromosome lengths from",infile)
    chromlengths=dict()
    for line in open(infile, "r"):
        data=line.strip().split()
        chrom , cum = data[0], data[2]
        if chrom.isnumeric():
            chrom = int(chrom)
        chromlengths[chrom] = int(cum)
    return chromlengths


def calc_empirical_pvals(results, permfile):
    perms = pd.read_table(permfile)

    # Convert to a dictionary of sorted numpy arrays for easier handling
    permvals = dict()
    for c in perms.columns:
        newname = re.sub(string=c, pattern="^chr", repl="")
        newname = int(newname)
        perm_pvals = np.array(sorted(perms[c]))  # Sort first to make p-val calc easier
        perm_pvals = perm_pvals[~np.isnan(perm_pvals)]  # Remove NAs
        permvals[newname] = perm_pvals

    # Go through and assign empirical p-values to each one
    for r in results:
        mychrom =np.unique(r['Chr'])
        if len(mychrom) != 1:
            print("WARNING!!! Trying to get p-values by chromosome but not 1 chromosome found:", mychrom)
        mychrom = mychrom[0]    # Get out of a list
        r['empirical_pval'] = np.searchsorted(a=permvals[mychrom], v = r['p'], side='left') / len(permvals[mychrom])

    return(results)


def plot_gwas(results, chrom_offsets, outprefix, rasterize=True, sparsify=None, trait="unknown", qtl=None, pval_cutoff=0.01):
    print("Plotting Manhatten plots with output prefix", outprefix)

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
    results = results.loc[~np.isnan(results['p']),]

    # Colors and size
    colors = 'black'
    size=5
    if 'empirical_pval' in results.columns:
        empiricals = np.array(results['empirical_pval'])   # To speed calculations up
        colors = np.array(['darkgray'] * len(results))
        colors[empiricals <= 0.1] = 'black'
        colors[empiricals <= 0.05] = 'blue'
        colors[empiricals <= 0.01] = 'darkred'

        # Adjust size
        size = np.array([5] * len(results))
        size[empiricals <= 0.1] = 15
        size[empiricals <= 0.05] = 20
        size[empiricals <= 0.01] = 30
    results['size'] = size
    results['color'] = colors

    # x-values
    xvals = list()
    for mychrom, mypos in zip(results['Chr'], results['Pos']):
        xvals.append(mypos + chrom_offsets[mychrom])
    results['x'] = xvals

    # Plot
    fig = plt.figure(figsize=(10, 3))
    grid = gridspec.GridSpec(nrows=100, ncols=100, hspace=0, wspace=0)
    ax_manhatten = fig.add_subplot(grid[5:95, :60], title=trait, xlabel="Position", ylabel="-log10 pvalue")
    myscatter = ax_manhatten.scatter(x=results['x'], y=-np.log10(results['p']), s=results['size'], color=results['color'],
                                     alpha=0.5, linewidths=0)
    myscatter.set_rasterized(rasterize)

    # Add chromosome divisions
    for c in chrom_offsets:
        ax_manhatten.axvline(chrom_offsets[c], color="lightgray", zorder=-1)

    # Add QTL locations
    if qtl is not None:
        qtl_pos = [mypos + chrom_offsets[mychrom] for mychrom, mypos in zip(qtl['chrom'], qtl['pos'])]
        for q in qtl_pos:
            ax_manhatten.axvline(q, color="purple", zorder=-2, linewidth=0.1)
        hits = results.loc[ results['empirical_pval'] <= pval_cutoff, :]
        for x in hits['x']:
            if x in qtl_pos: ax_manhatten.axvline(x, color="red", zorder=-1, linewidth=0.5)     # Any hits done in bright red


    # Add legend
    dot01 = patches.Patch(color='darkred', label='p ≤ 0.01')
    dot05 = patches.Patch(color='blue', label='p ≤ 0.05')
    dot10 = patches.Patch(color='black', label='p ≤ 0.10')
    ax_manhatten.legend(handles=[dot01, dot05, dot10], loc='upper right', fontsize='xx-small', frameon=False, title="Empirical\nP-values" )

    # Adjust figure size
    min_x, max_x = min(results['x']), max(results['x'])
    pad = (max_x - min_x)/10
    ax_manhatten.set_xlim(left=min_x-pad, right=max_x+pad*2)


    fig.savefig(outprefix + ".png", dpi=150)


def load_qtl(qtlfile, hapmap, trait):
    print("Loading QTL simulation parameters from",qtlfile)
    qtl_data = pd.read_table(qtlfile)
    qtl_data['name'] = [n.replace(".", "_") for n in qtl_data['name']]  # Change name to have just underscores so matche
    target = qtl_data['name'] == trait # Identify which row corresponds to this trait
    myqtl = qtl_data.loc[target, "qtn"].iloc[0] # Pull out comma-separated list of QTN names
    myqtl = set(myqtl.split(','))

    print("\tMatching with hapmap data from", hapmap)
    hmp = pd.read_table(hapmap)
    tokeep = [s in myqtl for s in hmp['rs#']]   # Identify which sites in the hapmap are actually tagged as QTL for this trait
    qtl = hmp.loc[tokeep, ['rs#', 'chrom', 'pos']]
    print("\t\tLocated",len(qtl),"qtl to plot")
    return qtl


def set_color(p, cutoff, chrom):
    if p < cutoff and chrom %2 == 0: return "blue"
    if p < cutoff and chrom %2 != 0: return "red"
    if chrom % 2 == 0: return "gray"
    if chrom % 2 != 0: return "silver"
    return "purple" # Should never get here

if __name__ == '__main__': main()