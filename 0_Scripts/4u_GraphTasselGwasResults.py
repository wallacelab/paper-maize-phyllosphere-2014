__author__ = 'jgwall'

import argparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import re

debug = False

matplotlib.rcParams.update({"font.size":'10'})

def main():
    args = parse_args()
    chromlengths = load_chromlengths(args.chromlengths)
    results = load_results(args.infiles)
    cutoffs = load_perms(args.permfile, args.perm_cutoff)
    plot_results(results, chromlengths, args.graphfile, args.rasterize, args.nmarkers, args.num_traits, cutoffs)
    output_top_hits(results, args.outfile, args.top_markers, cutoffs)


class gwas:
    name = None
    h2 = None
    chroms, poses, pvals = None, None, None

    def __init__(self, name, geneticvar, residvar):
        self.name = name
        self.h2 = geneticvar / (geneticvar + residvar)
        self.chroms, self.poses, self.pvals = list(), list(), list()

    def add_result(self, chr, pos, pval):
        self.chroms.append(chr)
        self.poses.append(pos)
        self.pvals.append(pval)

    def to_numpy(self):
        self.chroms = np.array(self.chroms)
        self.poses = np.array(self.poses)
        self.pvals = np.array(self.pvals)

    # Return the top n hits
    def top_n(self, n):
        order = np.argsort(self.pvals)
        best=order[:n]
        best_chroms = self.chroms[best]
        best_poses = self.poses[best]
        best_pvals = self.pvals[best]
        # Create a list of lines of text to output
        output=list()

        for chrom, pos, pval in zip(best_chroms, best_poses, best_pvals):
            output.append("\t".join([self.name, str(self.h2), str(chrom), str(pos), str(pval)]) + "\n")
        return output

    # Return hits below a certain cutoff
    def hits_below_cutoff(self, cutoffs):
        best, mycutoffs = list(), list()
        for i in range(len(self.pvals)):
            # print(self.pvals[i], "\t", self.chroms[i], "\t" , cutoffs[str(self.chroms[i])])
            if self.pvals[i] < cutoffs[str(self.chroms[i])]:
                best.append(i)
                mycutoffs.append(cutoffs[str(self.chroms[i])])
        best=np.array(best)

        if len(best) ==0: return ""

        best_chroms = self.chroms[best]
        best_poses = self.poses[best]
        best_pvals = self.pvals[best]
        # Create a list of lines of text to output
        output=list()

        for chrom, pos, pval, cut in zip(best_chroms, best_poses, best_pvals, mycutoffs):
            output.append("\t".join([self.name, str(self.h2), str(chrom), str(pos), str(pval), str(cut)]) + "\n")
        return output

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-g", "--graphfile")
    parser.add_argument("-t", "--top-markers", type=int, default=10, help="Write top N markers to output file")
    parser.add_argument("-c", "--chromlengths", help="File of chromosome lengths")
    parser.add_argument("-n", "--nmarkers", type=int, help="Number of markers in the original analysis (useful for making correct QQ plots)")
    parser.add_argument("-r", "--rasterize", default=False, action="store_true", help="Whether to rasterize the inner portion of each plot. (Only applicable to vector output formats)")
    parser.add_argument("-x", "--num-traits", type=int, help="Graph only the first X traits")
    parser.add_argument("-p", "--permfile", help="File of permuted GWAS results for doing null distribution")
    parser.add_argument("--perm-cutoff", type=float, default=0.05, help="Empirical error rate cutoff from permutations")
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



def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def load_results(infiles):

    n, skipped = 0, 0
    results=dict()
    for infile in infiles:
        print("Loading GWAS results from",infile)
        IN = open(infile, "r")

        # Parse header
        header =IN.readline().strip().split('\t')
        traitID, chrID, posID, pvalID, genvarID, residvarID = [header.index(s) for s in ["Trait", "Chr","Pos","p","Genetic Var","Residual Var"]]

        # print([traitID, chrID, posID, pvalID, genvarID, residvarID])

        # Go through and load results into a dictionary
        for line in IN:
            n+=1
            if n % 100000 == 0: print("\tProcessed",n,"lines")
            if debug and n >1000: break # For debugging
            data=line.strip().split('\t')
            trait = data[traitID]
            # Load trait if needed, including heritability value
            if trait not in results:
                results[trait]=gwas(trait, float(data[genvarID]), float(data[residvarID]))
            # Add result
            if data[chrID].isnumeric() and data[posID].isnumeric() and is_float(data[pvalID]):
                results[trait].add_result(chr=int(data[chrID]), pos=int(data[posID]), pval=float(data[pvalID]))
            else:
                skipped +=1
                n-=1

    # Solidify
    print("\tLoaded",len(results),"traits of data from",n+skipped,"lines of results (",skipped,"were skipped due to not being numeric. This should just be 1 per trait per input file)")
    for trait in results:
        print("\t\t",trait,"has",len(results[trait].poses),"datapoints on",len(np.unique(results[trait].chroms)),"different chromosomes")
        results[trait].to_numpy()
    return results


def plot_results(results, chromlengths, outfile, rasterize=False, nmarkers=None, numtraits=None, cutoffs=None):
    print("Plotting Manhatten plots to",outfile)
    # Set up plot parameters
    nplots = len(results) if numtraits is None else numtraits   # Set to user-specified # of traits if given, or just the total # of traits loaded
    width = 10
    height = 3 * nplots
    fig = plt.figure(figsize=(width, height))
    grid = gridspec.GridSpec(nrows=nplots, ncols=100, hspace=0.5, wspace=0.5)

    # Plot things
    row=0
    for trait in sorted(results.keys()):
        # Manhatten plot
        ax_manhatten = fig.add_subplot(grid[row,:60], title=trait, xlabel="Position", ylabel="-log10 pvalue")
        mytrait = results[trait]
        xvals = mytrait.poses + np.array([chromlengths[c] for c in mytrait.chroms])
        for mychrom in np.unique(mytrait.chroms):
            toplot = mytrait.chroms == mychrom
            pvals = mytrait.pvals[toplot]
            pvals[np.isnan(pvals)]=1
            color = "blue" if mychrom % 2 == 0 else "red"
            if cutoffs:
                color = [set_color(p, cutoffs[trait][str(mychrom)], mychrom) for p in pvals]
            myscatter = ax_manhatten.scatter(xvals[toplot], -np.log10(pvals), s=1, color=color)
            myscatter.set_rasterized(rasterize)

        for c in chromlengths:
            ax_manhatten.axvline(chromlengths[c], color="lightgray", zorder=-1)

        # QQ plot
        ax_qq = fig.add_subplot(grid[row,75:], title="QQ plot", xlabel="null distribution", ylabel="Actual distribution")
        if nmarkers is None:
            nmarkers = len(xvals)
        actual = np.sort(mytrait.pvals)
        null = np.arange(1, nmarkers+1) / nmarkers
        # Cut null down so matches actual values (assumed to be the most significant)
        if len(actual) < len(null):
            null = null[:len(actual)]
        # Reverse values so plot right
        actual = actual[::-1]
        null = null[::-1]
        myqq = ax_qq.scatter(-np.log10(null), -np.log10(actual), s=2, color="black")
        myqq.set_rasterized(rasterize)
        # Draw line
        xlim, ylim = ax_qq.get_xlim(), ax_qq.get_ylim()
        ax_qq.plot([0,10], [0,10], 'b-')
        ax_qq.set_xlim(left=xlim[0], right=xlim[1]) # so doesn't resize from the plotted line
        ax_qq.set_ylim(bottom=ylim[0], top=ylim[1])
        

        row+=1
        if row >= nplots: break
    fig.savefig(outfile, dpi=600)

def output_top_hits(results, outfile, n, cutoffs):
    print("Writing top",n,"results for each trait to",outfile)
    OUT = open(outfile, "w")
    header = ["trait","h2","chrom","pos","pval"]
    if cutoffs: header+= ['p_cutoff']
    OUT.write("\t".join(header) + "\n")    # Header
    for trait in sorted(results.keys()):
        if cutoffs is None:
            OUT.writelines(results[trait].top_n(n))
        else:
            OUT.writelines(results[trait].hits_below_cutoff(cutoffs[trait]))
    OUT.close()

    # Write cutoffs out
    if cutoffs:
        outfile = re.sub(string=outfile, pattern=".txt", repl=".cutoffs.txt")
        print("Writing p-value cutoffs to",outfile)
        OUT = open(outfile, "w")
        OUT.write("chrom\tcutoff\ttrait\n")
        for trait in sorted(results.keys()):
            for chrom in sorted(cutoffs[trait].keys()):
                OUT.write("\t".join([str(chrom), str(cutoffs[trait][chrom]), str(trait)]) + '\n')
        OUT.close()


def load_perms(permfile, cutoff):
    if permfile is None:
        print("No permutation file found; raw results reported")
        return None
    print("Loading permuted GWAS results from",permfile,"with cutoff",cutoff)

    permdata = pd.read_csv(permfile, sep='\t')
    perms = dict()

    # Split by trait
    traits = set(permdata["trait"])
    for t in traits:
        subset = permdata.loc[permdata['trait']==t,:]
        perms[t]=dict()
        for chrom in subset.columns:
            if chrom == "trait" : continue
            perms[t][chrom] = np.array(subset[chrom])
    #     print("\tLoaded",len(perms[t]),"chromosomes for",t)
    # print("\tLoaded", len(perms), "traits of permutations")

    # Calculate cutoffs based on permuted p-values on a per-chromosome basis
    cutoffs=dict()
    for t in perms:
        cutoffs[t]=dict()
        for chrom in perms[t]:
            pvals = sorted(perms[t][chrom])
            cutoff_index = int(cutoff * len(pvals))
            if cutoff_index >= len(pvals): cutoff_index = len(pvals) -1 # to prevent out-of-bounds
            cutoffs[t][chrom] = pvals[cutoff_index] # Cutoff value; p-values _less than_ this (not "<=" ) are significant

    return cutoffs

# Helper function to set color based on p-value and chromosome
def set_color(p, cutoff, chrom):
    if p < cutoff and chrom %2 == 0: return "blue"
    if p < cutoff and chrom %2 != 0: return "red"
    if chrom % 2 == 0: return "gray"
    if chrom % 2 != 0: return "silver"
    return "purple" # Should never get here

if __name__ == '__main__': main()