__author__ = 'jgwall'

import argparse

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

debug = False


def main():
    args = parse_args()
    print("Plotting combined GWAS hits for clustered SNPs in",args.infile)

    # Load data
    chrom_offsets = load_chromlengths(args.chromlengths)
    data = pd.read_table(args.infile)
    data = data.sort(['Chr', 'Pos', 'Trait'])   # Sort
    data[['Chr', 'Pos']] = data[['Chr', 'Pos']].astype(int, errors="raise") # convert floats to ints
    print("\tLoaded",len(data),"sites across", len(np.unique(data['Trait'])),"traits")

    # Break into chromosome windows for plotting
    plotdata = break_into_windows(data, args.winsize)

    # Output text
    outtext = args.outprefix + ".txt"
    output = pd.concat(plotdata)
    output['group_pos'] = output['group_pos'].astype(int)
    print("Outputting combined text file to", outtext)
    output.to_csv(outtext, sep='\t', index=None)

    # # Plot
    cors = pd.read_table(args.cors)
    plot_gwas(plotdata, chrom_offsets, args.outprefix, args.cutoffs, args.winsize, correct_for_correlation=False)
    plot_gwas(plotdata, chrom_offsets, args.outprefix, args.cutoffs, args.winsize, correct_for_correlation=True, cors=cors)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Input file of compiled hits for all GWAS traits")
    parser.add_argument("-o", "--outprefix", help="Output file prefix")
    parser.add_argument("-c", "--cutoffs", nargs="*", type=float, default=[0.1, 0.05, 0.01], help="List of cutoffs for graphing empirical p-values")
    parser.add_argument("-w", "--winsize", type=int, default=100000, help="How big to make the window size for plotting")
    parser.add_argument("--chromlengths", help="File with chromosome lengths and their offsets")
    parser.add_argument("--cors", help="File with matrix of correlations among all the BLUPs (for correcting the counts in each window")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def load_chromlengths(infile):
    print("Loading chromosome lengths from",infile)
    chromlengths=dict()
    mymax=-1    # Maximum value
    for line in open(infile, "r"):
        data=line.strip().split()
        chrom , length, cum = data
        if chrom.isnumeric(): chrom = int(chrom)
        chromlengths[chrom] = int(cum)
        mymax = max(mymax, int(cum) + int(length))
    chromlengths["MAX"] = mymax
    return chromlengths


def break_into_windows(data, winsize):
    print("Breaking hits into windows of", winsize)

    # First, break by chromosome
    data['group_pos'] = -1
    chroms = data.groupby('Chr')

    # Now iterate over and break into windows within each chromosome
    splitdata = list()
    for mychrom, mydata in chroms:
        mydata = mydata.copy()  # To avoid warnings of this being a subslice of a larger dataframe
        min = 0
        max = np.max(mydata['Pos'])
        splits = np.arange(start=min, stop=max + winsize, step=winsize)
        print("\tChrom",mychrom,"has",len(splits),"windows between min",min,"and max",max)


        # There's probably a python library to do this next part of identifying which split each goes in,
        # but the cost of learning is higher than just brtue-forcing it
        for i in range(len(splits)-1):
            group_pos = (splits[i] + splits[i+1]) / 2
            targets = (mydata['Pos'] >= splits[i]) & (mydata['Pos'] < splits[i+1])
            if np.sum(targets) == 0 : continue
            if np.any(mydata['group_pos'].loc[targets] != -1):
                print("\tWARNING! Group starting at",splits[i],"contains data points that have already been assigned a group")
            mydata.loc[targets, 'group_pos'] = group_pos

        #Split by window within chromosome
        subsplits = mydata.groupby("group_pos")
        for subname, subgroup in subsplits:
            splitdata.append(subgroup.copy())
    print("\tTotal",len(splitdata),"informative windows identified")
    return splitdata


def plot_gwas(data, chrom_offsets, outprefix, cutoffs, winsize, cors=None, correct_for_correlation=False):
    print("Plotting compiled GWAS results for empirical p-value cutoffs of",cutoffs)

    # Set up figure
    ncol=1
    nrow=len(cutoffs)
    fig = plt.figure(figsize=(ncol * 8, nrow * 4))
    grid = gridspec.GridSpec(ncols=ncol, nrows=nrow)
    hitdata = list()

    # Plot each cutoff level as a new axes
    index_i=0   #for compiling hit data as a DataFrame
    for row in range(nrow):
        cutoff = cutoffs[row]
        ax = fig.add_subplot(grid[row, :], title="GWAS results binned at " + str(winsize) + " and empirical pval cutoff "+str(cutoff))

        # Determine count in each window
        xvals, yvals, = list(), list()
        for window in data:

            # Calculate y values
            count, unique_traits = get_gwas_count(window, cutoff, cors, correct_for_correlation)
            if count == 0: continue
            yvals.append(count)

            # Calculate x values
            mychrom = np.unique(window['Chr'])
            mypos = np.unique(window['group_pos'])
            if len(mychrom) > 1: print("WARNING ! More than one chromosome found in window", window)
            if len(mypos) > 1: print("WARNING ! More than one plot location found in window", window)
            mychrom, mypos = mychrom[0], mypos[0]
            xval = chrom_offsets[mychrom] + mypos
            xvals.append(xval)

            # Store counts of hits for later output
            index_i+=1
            mytraits = ",".join(unique_traits)
            tempdata = pd.DataFrame({"cutoff":cutoff, "chrom":mychrom, "window_pos":mypos, "num_traits":len(unique_traits), "corrected_count":count, "traits":mytraits}, index=[index_i])
            hitdata.append(tempdata)

        # Plot as ball-and-stick model
        ax.scatter(xvals, yvals)
        ax.vlines(xvals, ymin=0, ymax=yvals)

        # Prettify axis
        xpad = 1e7
        ax.set_xlim(left = -xpad, right = chrom_offsets['MAX'] + xpad)
        add_chrom_divisions(ax, chrom_offsets)
        # break

    # Save hit data
    outhits = pd.concat(hitdata)
    outhits = outhits[['cutoff','chrom','window_pos','num_traits','corrected_count', 'traits']]
    outhits = outhits.sort(['cutoff','chrom','window_pos'])
    outhits.to_csv(outprefix+".hitcounts.txt", sep='\t', index=None)

    # Save image
    outfile=outprefix + ".png"
    if correct_for_correlation: outfile = outprefix + ".corrected.png"
    print("\tSaving figure to", outfile)
    fig.savefig(outfile, dpi=100)


# Helper function to determine how many hits to return
def get_gwas_count(window, cutoff, cors, correct_for_correlation=False):
    mydata = window.loc[window['empirical_pval'] <= cutoff,:]
    unique_traits = sorted(np.unique(mydata['Trait']))
    count = len(unique_traits)  # Default; each trait adds a value of 1

    # If need to correct for correlation, do so so that each count is actually (1- highest correlation among earlier traits)
    if correct_for_correlation and len(unique_traits) > 0:
        mycounts = [0] * len(unique_traits)
        for i in range(len(unique_traits)):
            # If first trait, assign a value of 1 and continue
            if i==0:
                mycounts[i]=1
                continue

            mytrait=unique_traits[i]
            prev_traits = unique_traits[:i]
            max_cor = find_max_correlation(mytrait, prev_traits, cors)
            mycounts[i] = 1-max_cor

            # print("Mytrait:", mytrait)
            # print("Previous traits:", prev_traits)
            # print('\n')

        # print("#####\n",mycounts)
        count = sum(mycounts)
        # print(count,"\n#####\n")

    # # For debugging
    # print(window)
    # print("Find",count,"at cutoff",cutoff,"\n")

    return count, unique_traits

def add_chrom_divisions(ax, offsets):
    for o in offsets:
        ax.axvline(offsets[o], linestyle='dashed', color='gray', zorder=-99)

# helper function to find the maximum correlation between a trait and other listed traits
def find_max_correlation(target, others, cors):
    subset = cors.loc[target, others]
    mycors = [abs(c) for c in subset]
    mymax=max(mycors)
    # print(subset)
    # print(mycors, mymax)
    return(mymax)

if __name__ == '__main__': main()