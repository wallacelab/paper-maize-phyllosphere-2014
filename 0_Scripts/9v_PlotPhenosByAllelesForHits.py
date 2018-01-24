__author__ = 'jgwall'

import argparse
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

debug = False


def main():
    args = parse_args()
    print("Plotting phenotypes by alleles for hits in",args.infile,"at p-value cutoff of",args.pval_cutoff)
    hits=pd.read_table(args.infile)
    genos=pd.read_table(args.genos)
    phenos=pd.read_table(args.phenos, skiprows=2, index_col="Taxa")

    # Subset based on empirical p-value
    hits=hits.loc[hits['empirical_pval'] <= args.pval_cutoff]
    hits = hits.sort(columns=['Trait','Chr','Pos'])

    # Match up genotypes and phenotypes
    alleles = genos.iloc[:,11:]
    alleles.index = genos['rs#']
    phenos = phenos.loc[alleles.columns,:]
    if len(phenos) != len(alleles.columns):
        print("\tWARNING! Genotype and phenotype are of different lengths after matching")
    mismatch = [myp != myg for myp, myg in zip(phenos.index, alleles.columns)]
    if np.any(mismatch):
        print("\tWARNING!",np.sum(mismatch),"mismatches found after trying to match genotype and phenotype")

    # Determine number of plots
    nplots = len(hits)
    traits, traitcount = np.unique(hits['Trait'], return_counts=True)
    nrow = len(traits)
    ncol = max(traitcount)
    print("Plotting",nplots,"hits in a ",nrow,"x",ncol,"grid")

    # Set up figure
    fig = plt.figure(figsize=[ncol*3, nrow*3])
    grid =gridspec.GridSpec(nrows=nrow, ncols=ncol, hspace=0.5, wspace=0.5)

    # Loop over and plot
    myrow, mycol = 0, 0
    for i in range(len(hits)):
        mytrait = hits['Trait'].iloc[i]
        mymarker = hits['Marker'].iloc[i]

        # Plot violin plots of data
        ax = fig.add_subplot(grid[myrow, mycol])
        labels, values = get_plotdata(mytrait, mymarker, phenos, alleles)
        xticks=np.arange(len(labels))
        ax.violinplot(dataset=values, positions=xticks)

        # Label
        ax.set_title(mytrait + "\n" + mymarker, weight='bold', fontsize='x-small')
        ax.set_xticks(xticks)
        ax.set_xticklabels(labels)
        ax.set_xlim(left=min(xticks)-1, right=max(xticks)+1)

        # Add number of observations below each violin
        for group in xticks:
            ax.text(x=group, y=ax.get_ylim()[0], s=len(values[group]), weight='bold', fontsize='xx-small', ha='center', va='bottom')


        # Increment row/col
        mycol +=1
        if i < len(hits)-1 and hits['Trait'].iloc[i+1] != mytrait:  # Advance row when hit a new trait
            mycol=0
            myrow +=1


    fig.savefig(args.outfile, dpi=100)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="File of clustered SNP hits from 4x_ClusterSnpHits.r (with the .best.txt suffix)")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-g", "--genos", help="Hapmap of SNPs that will be needed")
    parser.add_argument("-P", "--phenos", help="TASSEL-formatted phenotype file")
    parser.add_argument("-p", "--pval-cutoff", type=float, default=1, help="Empirical p-value cutoff for including hits")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def get_plotdata(trait, marker, phenos, alleles):
    myphenos = np.array(phenos.loc[:,trait])
    mygenos = np.array(alleles.loc[marker, :])

    plotdata=dict()
    for pheno, allele in zip(myphenos, mygenos):
        if allele == "N" or allele == "NN": continue    # Skip missing data
        if np.isnan(pheno): continue    # Skip missing phenotype data
        if allele not in plotdata:
            plotdata[allele] = list()
        plotdata[allele].append(pheno)

    # Get lists of vlaues and labels for the plot (since violinplot() doesn't take a dict
    values, labels = list(), list()
    for allele in sorted(plotdata.keys()):
        labels.append(allele)
        values.append(plotdata[allele])
    return labels, values


if __name__ == '__main__': main()