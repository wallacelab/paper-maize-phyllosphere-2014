__author__ = 'jgwall'

import argparse
import itertools
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind as t_test
import string

debug = False


def main():
    args = parse_args()
    print("Calculating day-night distances from",args.infile)
    dist, labels = load_distances(args.infile)
    np.random.seed(args.seed)

    # Get pairs of samples to test
    pairs = get_sample_pairs(dist, labels)  # Sample plot day-night pairs
    neighbors = get_neighbors(dist, labels) # Nieghboring plots, day-night pairs
    randoms = get_randoms(dist, labels, pairs, neighbors)      # Random plots, day-night pairs

    # Plot
    fig = plt.figure(figsize=(10,5))
    grid = gridspec.GridSpec(nrows=10, ncols=3, wspace=0.2)
    ylim = [0, 1.05] if ".weighted_unifrac" not in args.outprefix else [0, 1.25]    # Weighted unifract distances can be >1, so have to adjust
    ax_8 = fig.add_subplot(grid[:8,0], title="August 8 Samples")
    ax_26 = fig.add_subplot(grid[:8, 1], title="August 26 Samples")
    ax_all = fig.add_subplot(grid[:8, 2], title="All Samples")
    plot_violins(ax_8, pairs, neighbors, randoms, date="8", pval_outfile=args.outprefix+".pvals_aug8.txt", ylim=ylim)
    plot_violins(ax_26, pairs, neighbors, randoms, date="26", pval_outfile=args.outprefix+".pvals_aug26.txt", ylim=ylim)
    plot_violins(ax_all, pairs, neighbors, randoms, date=None, pval_outfile=args.outprefix+".pvals_all.txt", ylim=ylim)


    fig.savefig(args.outprefix +".png", dpi=100)
    fig.savefig(args.outprefix + ".svg", dpi=100)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="QIIME distance matrix file")
    parser.add_argument("-o", "--outprefix", help="Output file prefix")
    parser.add_argument("-s", "--seed", type=int, default=1, help="Random seed")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def load_distances(infile):
    dist = pd.read_table(infile, index_col=0)
    rownames, colnames = dist.index, dist.columns
    mismatches = [r!=c for r,c in zip(rownames, colnames)]
    if any(mismatches):
        print("WARNING! Row and column labels do not match!")

    labels=dict()
    for i in range(len(rownames)):
        sample = rownames[i]
        labels[sample]=i

    dist = np.array(dist)
    return dist, labels

def get_sample_pairs(dist, labels):
    print("\tExtracting day-night pairs from the same samples")
    pairs=dict()
    n=0
    for s1 in sorted(labels.keys()):
        s2=switch_day_night(s1)
        if s2 not in labels: continue
        pairs = add_to_sampleset(pairs, s1, s2, labels, dist)
        n+=1
    print("\t\tIdentified",n,"total pairs that collapse to",len(pairs),"unique pairs of the same plot in day-night combinations")
    return pairs

def get_neighbors(dist, labels):
    print("\tExtracting day-night pairs from neighboring samples")
    pairs=dict()
    n=0
    for s1 in sorted(labels.keys()):
        s2 = switch_day_night(s1) # Change from day to night and vice-versa
        # Get neighboring by adding or subtracting 2 from the plot number
        s2_left = change_plot(s2, -2)
        s2_right = change_plot(s2, 2)

        # Check if corresponding samples exist, and add
        if s2_left in labels:
            pairs = add_to_sampleset(pairs, s1, s2_left, labels, dist)
            n+=1
        if s2_right in labels:
            pairs = add_to_sampleset(pairs, s1, s2_right, labels, dist)
            n+=1
    print("\t\tIdentified",n,"total pairs that collapse to",len(pairs),"unique pairs of neighboring plots in day-night combinations")
    return pairs

def get_randoms(dist, labels, sample_pairs, neighbors):
    print("\tExtracting day-night pairs from random samples that are not in the above sets")
    pairs=dict()
    n=0
    samples = sorted(labels.keys())
    for i in range(len(samples)):
        s1 = samples[i]
        for j in range(i, len(samples)):
            s2=samples[j]

            if s1 == s2: continue   # If samples are same, skip. (Shouldn't happen, but best to be safe)

            # If both are day or both are night, skip
            if s1.startswith("LMAN") and s2.startswith("LMAN"): continue
            if s1.startswith("LMAD") and s2.startswith("LMAD"): continue

            # If sample pair is part of either previous set, skip
            key = make_key(s1, s2)
            if key in sample_pairs or key in neighbors: continue

            # If samples are on different dates, skip
            if get_date(s1) != get_date(s2):
                # print("Plots",s1,s2,"are on different dates and so skipping")
                continue
            # else: print("Plots",s1,s2,"are on the same date")

            # Add pairs that made it this far to the set
            pairs = add_to_sampleset(pairs, s1, s2, labels, dist)
            n+=1
    print("\t\tIdentified",n,"total pairs that collapse to",len(pairs),"unique pairs of random plots in day-night combinations on same date")
    return pairs


# Do string substitution to get corresponding day or night value
def switch_day_night(s1):
    if "LMAN" in s1: s2 = s1.replace("LMAN", "LMAD")
    if "LMAD" in s1: s2 = s1.replace("LMAD", "LMAN")
    return s2

def change_plot(sample, increment):
    plot = int(sample[-4:])
    # print("Sample",sample,"has plot number",plot)
    plot += increment
    newsample = sample[:-4] + "{0:0>4}".format(plot)    # formatting is just how to get leading zeroes okay
    # print("\tNeighbor has sample ID",newsample)
    return newsample

# Extract distance and add to a dictionary of
def add_to_sampleset(pairs, s1, s2, labels, dist):
    row, col = labels[s1], labels[s2]
    mydist = dist[row, col]
    # print("\t\tPair", s1, s2, "is at", row, col, "with distnace", mydist)

    # Make unique by sorting and adding to dictionary
    key = make_key(s1, s2)
    if key in pairs:    # If already added, check that distance values are the same
        if pairs[key] != mydist:
            print("\t\tWARNING!! Pair",key,"already added but previous distance",pairs[key],"!= new distance",mydist)
        else: pass
    else:   # If not present, add
        pairs[key] = mydist

    return pairs

def make_key(s1, s2):
    return "-".join(sorted([s1, s2]))

def unmake_key(key):
    return key.split("-")

def get_date(sample):
    return sample.split('.')[1]


def plot_violins(ax, pairs, neighbors, randoms, date=None, pval_outfile=None, ylim=[0, 1.05]):
    print("Plotting violins for date",date)
    if date is not None:    # Subset by date if required
        pairs = subset_pairs_by_date(pairs, date)
        neighbors = subset_pairs_by_date(neighbors, date)
        randoms = subset_pairs_by_date(randoms, date)

    # Randomly subset the randoms to same size as the pair data
    random_subset = np.random.choice(sorted(randoms.keys()), size=len(pairs))
    subset_data = [randoms[k] for k in random_subset]

    # Plot
    pair_data, neighbors_data, randoms_data = list(pairs.values()), list(neighbors.values()), list(randoms.values())
    plot_data = [pair_data, neighbors_data, randoms_data, subset_data]
    xvals = [1,2,3,4]
    violins = ax.violinplot(plot_data, positions=xvals)

    # Statistical tests
    lettercodes, pvals = get_significant_groupings(plot_data, labels=["Pairs","Neighbors","Randoms","Randoms_subset"], cutoff=0.01)
    maxes = [max(d) for d in plot_data]
    for x,y,s in zip(xvals, maxes, lettercodes):
        ax.text(x=x, y=y+0.02, s=s, horizontalalignment='center', fontweight='bold', fontsize='small')
    print("P-values for significance tests:\n",pvals)

    # Write out if requested
    if pval_outfile:
        pvals['dataset'] = ax.get_title()
        pvals.to_csv(pval_outfile, sep='\t')


    # Clean up axes
    ax.set_xticks(xvals)
    ax.set_xticklabels(["Same plots","Neighbor plots","Random plots","Random plots\n(subset)"], rotation="vertical", fontsize="x-small", fontweight="bold")
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('left')
    ax.set_title(ax.get_title(), fontsize="large", fontweight="bold")

    # Y axis stuff
    ax.set_ylim(ylim)
    yticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    yticklabels=[str(y) for y in yticks]
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize='small', weight='bold')

    # Clean up violins
    for body in violins['bodies']:
        body.set_color('cornflowerblue')
        body.set_alpha(1)
        body.set_linewidth(0)
    for component in ['cmins','cmaxes','cbars']:
        violins[component].set_color('midnightblue')

def subset_pairs_by_date(pairs, date):
    newpairs=dict()
    for key in pairs:
        s1, s2 = unmake_key(key)
        d1, d2 = get_date(s1), get_date(s2)
        if d1 != d2: print("\t\tWARNING! Dates are different for pair",key)

        # If sample pair was on the indicated date, keep. Otherwise discarded
        if d1 == date:
            newpairs[key] = pairs[key]
    print("\t",len(newpairs),"of",len(pairs),"sample pairs kept for date",date)
    return(newpairs)

# function to take a list of datapoints in a group, do pairwise t-tests, determine which groups are statistically significant, and return the letter codes for groupings
def get_significant_groupings(groupdata, labels=None, cutoff=0.05):
    print("\t\tDetermining significant letter groupings for plot\n")
    if labels==None: labels=["Data" + str(i+1) for i in range(len(groupdata))]   # Labels for each dataset
    
    # Calculate p-values
    pvals=np.empty(shape=(len(groupdata), len(groupdata)))
    pvals[:] = np.nan
    for i in range(len(groupdata)):
        myset = set()
        for j in range(i+1, len(groupdata)):
            t, p = t_test(groupdata[i], groupdata[j])
            pvals[i,j], pvals[j,i] = p, p
    pvals=pd.DataFrame(pvals, columns=labels, index=labels)
    print("P-values:\n",pvals)

    # Make all possible combinations of groups
    allgroups=list()
    for size in range(2, len(labels)+1):
        for mygroup in itertools.combinations(labels, size):
            allgroups.append(mygroup)
            #print(mygroup)
    
    # Subset to just the ones where all members are statistically indistinguishable
    goodgroups=list()
    for mygroup in allgroups:
        isgood=True
        for i in range(len(mygroup)):
            for j in range(i, len(mygroup)):
                if pvals.loc[mygroup[i], mygroup[j]] <= cutoff:
                    isgood=False
        #print(isgood, mygroup)
        if isgood: goodgroups.append(mygroup)
    
    # Check if one-member groups exist by collapsing groups to a set of labels
    already_included = set(np.hstack(goodgroups))
    for l in labels:
        if l not in already_included: 
            goodgroups.append([l])  #Append 1-item list of groups
    #print(goodgroups)

    # Now assign letter values; do in a loop over labels to preserve order
    lettercodes = {l:"" for l in labels}    # Dictionary of which letter codes each label has
    letter_i=0  # Keeping track of which group label we're using
    for l in labels:
        for g in goodgroups:
            if l not in g: continue # Skip groups this label is not in
            for mylabel in g:
                lettercodes[mylabel] += string.ascii_lowercase[letter_i]
            goodgroups.remove(g)    # Remove a group after it's been processed so it doesn't get processed a second time
        letter_i+=1
    #print(lettercodes)
    
    # Make letters in order of original list
    lettervals = [lettercodes[l] for l in labels]
    #print(lettervals)
    return lettervals, pvals
    
    

if __name__ == '__main__': main()