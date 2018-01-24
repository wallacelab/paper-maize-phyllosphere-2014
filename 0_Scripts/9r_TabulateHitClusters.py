__author__ = 'jgwall'

import argparse

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

debug = False


def main():
    args = parse_args()
    print("Tabulating clusters of repeated hits in",args.infile)

    # Load data
    data = pd.read_table(args.infile)
    cors = pd.read_table(args.cors)
    data = data.sort(['Chr', 'Pos', 'Trait'])   # Sort
    data[['Chr', 'Pos']] = data[['Chr', 'Pos']].astype(int, errors="raise") # convert floats to ints
    print("\tLoaded",len(data),"sites across", len(np.unique(data['Trait'])),"traits")

    # Remove hits with too low of a p-value
    data = data.loc[data['empirical_pval'] <= args.pval_cutoff]
    print("\t", len(data), "remain after filtering for empirical p-values <=",args.pval_cutoff)

    # Break into chromosome windows for plotting
    mydata = break_into_windows(data, args.winsize)

    # Output raw results of which hits go in which region
    print("Outputting raw text file")
    parsed = pd.concat(mydata)
    parsed['window'] = ["chr" + str(c) + "_" + str(int(p)) for c, p in zip(parsed['Chr'], parsed['window'])]
    parsed = parsed.sort(columns=['Trait','Chr','Pos'])
    parsed.to_csv(args.outprefix + ".raw.txt", sep='\t', index=None)

    # Get report for each window
    windows = np.unique(parsed['window'])
    tmp_reports = list()
    traits, counts, scores = list(), list(), list()
    for mywin in windows:
        myhits = parsed.loc[parsed['window']==mywin].copy() # Subset to window of interest
        myscore, mytraits = get_gwas_count(myhits, cors)

        myhits['n_traits'] = len(mytraits)
        myhits['score'] = myscore
        tmp_reports.append(myhits)

    # Write out full report so can look at individual hits
    print("Outputting window reports")
    full = pd.concat(tmp_reports)
    full=full[['window','n_traits','score','Trait','Marker','Chr','Pos','p','empirical_pval']]
    full.to_csv(args.outprefix + ".windows_full.txt", sep='\t', index=None)

    # Collapse to abbreviated report
    windows, markers, counts, scores, traits, locs = list(), list(), list(), list(), list(), list()
    for mywin in full.groupby('window'):
        windows.append(mywin[0])    # Store value used to split
        mywin=mywin[1]  # now just work with the subsetted dataframe
        counts.append(np.mean(mywin['n_traits']))
        scores.append(np.mean(mywin['score']))
        traits.append(",".join(mywin['Trait']))
        markers.append(",".join(np.unique(mywin['Marker'])))

        # Locations
        mylocs = [str(c) + ":" + str(p) for c,p in zip(mywin['Chr'],mywin['Pos'])]
        mylocs = ",".join(np.unique(mylocs))
        locs.append(mylocs)

    short = pd.DataFrame({'window':windows, 'n_traits':counts,'score':scores,'locations':locs,'traits':traits, 'markers':markers})
    short=short[['window','n_traits','score','locations','markers', 'traits']]
    short.to_csv(args.outprefix + ".windows_short.txt", sep='\t', index=None)

    # Filter for just the high-quality clusters
    full = full.loc[full['score'] >= args.min_score,]
    short = short.loc[short['score'] >= args.min_score,]
    full.to_csv(args.outprefix + ".windows_full.filtered.txt", sep='\t', index=None)
    short.to_csv(args.outprefix + ".windows_short.filtered.txt", sep='\t', index=None)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Input file of compiled hits for all GWAS traits")
    parser.add_argument("-o", "--outprefix", help="Output file prefix")
    parser.add_argument("-w", "--winsize", type=int, default=100000, help="How big to make the window size for plotting")
    parser.add_argument("-p", "--pval-cutoff", type=float, default=0.05, help="Empirical p-value cutoff for including hits")
    parser.add_argument("-m", "--min-score", type=float, default=1, help="Minimum score value for going into the filtered output files")
    parser.add_argument("--cors", help="File with matrix of correlations among all the BLUPs (for correcting the counts in each window")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def break_into_windows(data, winsize):
    print("Breaking hits into windows of", winsize)

    # First, break by chromosome
    data['window'] = -1
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
            window = (splits[i] + splits[i+1]) / 2
            targets = (mydata['Pos'] >= splits[i]) & (mydata['Pos'] < splits[i+1])
            if np.sum(targets) == 0 : continue
            if np.any(mydata['window'].loc[targets] != -1):
                print("\tWARNING! Group starting at",splits[i],"contains data points that have already been assigned a group")
            mydata.loc[targets, 'window'] = window

        #Split by window within chromosome
        subsplits = mydata.groupby("window")
        for subname, subgroup in subsplits:
            splitdata.append(subgroup.copy())
    print("\tTotal",len(splitdata),"informative windows identified")
    return splitdata



# Helper function to determine how many hits to return
def get_gwas_count(window, cors):
    unique_traits = sorted(np.unique(window['Trait']))

    # Correct for correlation so that score is actually (1- highest correlation among earlier traits)
    mycounts = [0] * len(unique_traits)
    for i in range(len(unique_traits)):
        # If first trait, assign a value of 1 and continue
        if i==0:
            mycounts[i]=1
            continue
        # Otherwise correct for corrleations with previously added traits
        mytrait=unique_traits[i]
        prev_traits = unique_traits[:i]
        max_cor = find_max_correlation(mytrait, prev_traits, cors)
        mycounts[i] = 1-max_cor
        # print("Mytrait:", mytrait)
        # print("Previous traits:", prev_traits)
        # print('\n')
    score = sum(mycounts)
    return score, unique_traits

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