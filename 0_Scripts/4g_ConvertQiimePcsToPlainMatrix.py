__author__ = 'jgwall'

import argparse
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
from os.path import commonprefix
import re

debug = False

matplotlib.rcParams.update({'font.size':10})

def main():
    args = parse_args()
    pcs, eigens = load_pcs(args.infiles, args.num_pcs)
    output_matrix(pcs, args.outfile)
    if args.outgraphic:
        output_pheno_distributions(pcs, args.outgraphic, args.num_pcs)
    if args.screeplot:
        output_scree_plots(eigens, args.infiles, args.screeplot, args.num_eigens)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*", help="Principal components files from QIIME")
    parser.add_argument("-o", "--outfile", help="Output phenotype file in TASSEL format")
    parser.add_argument("-g", "--outgraphic", help="Output graphic of histrograms of each phenotype")
    parser.add_argument("-s", "--screeplot", help="Output graphic of scree plots for the PCs")
    parser.add_argument("-n", "--num-pcs", type=int, default=1, help="Number of principal components to take")
    parser.add_argument("--num-eigens", type=int, default=20, help="Number of eigenvalues to plot in the screeplot")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def load_pcs(infiles, num_pcs):
    if debug: infiles=infiles[:5]   # Debug mode = take only first 5 input files
    print("Loading principal components from",len(infiles),"input files")
    data = [load_pc_file(f) for f in infiles]
    eigens = [load_eigens(f) for f in infiles]
    print("\tTaking the first",num_pcs,"components of each")
    data = [d[:,:num_pcs+1] for d in data]

    # Print variance explained
    for eigenvals, file in zip(eigens, infiles):
        print("\tVariance explained by each PC:")
        for i in range(num_pcs):
            print("\t\tPC",i+1, ":",eigenvals[i])
        print("\t\tTOTAL:", sum(eigenvals[:num_pcs]))

    # Handle dataset names
    names = [i.replace('/', '_') for i in infiles]
    left_trim = commonprefix(names)
    right_trim = reverse(commonprefix([reverse(n) for n in names])) # Same as above, but doing it backwards


    newnames = [re.sub(pattern="^" + left_trim, repl="", string=s) for s in names]
    newnames = [re.sub(pattern=right_trim + "$", repl="", string=s) for s in newnames]
    # for old, new in zip(names, newnames):
    #     print(new,"from",old)

    # Convert to individual DataFrames for easier unification
    colnames = ["PC" + str(p) for p in range(1, num_pcs+1)]
    data = [pd.DataFrame(d[:,1:], index=d[:,0], dtype=float, columns=colnames) for d in data]  # Slice first columns out to be the index
    for mydata, myname in zip(data, newnames):
        mydata.columns = [myname + "_" + c for c in mydata.columns] # Add dataset name as part of column name
    # print(data[0].head())

    # Unify all DataFrames into one
    master = pd.concat(data, axis=1)
    # print(master.head())

    # Return a concatenated dataframe
    return master, eigens

# Helper function to reverse an iterable
def reverse(s):
    return s[::-1]

# Helper function to load a single PC file from QIIME
def load_pc_file(infile):
    IN=open(infile, "r")
    for line in IN:    # Advance through lines until get to right spot
        if line.startswith("Site"): break
    site, nrow, ncol = line.strip().split('\t')
    matrix = list()
    for line in IN:
        if line == "\n": break  # When hit white space, break out
        data = line.strip().split('\t')
        matrix.append(data)
    matrix = np.array(matrix)
    IN.close()
    return matrix

def load_eigens(infile):
    IN=open(infile, "r")
    eigens=None
    for line in IN:    # Advance through lines until get to right spot
        if line.startswith("Proportion explained"):
            eigens = IN.readline().strip().split('\t')
            eigens = [float(e) for e in eigens]
            break
    IN.close()
    return eigens

def output_matrix(pcs, outfile):
    print("Outputting",len(pcs.columns), "PCs as matrix to",outfile)
    pcs.to_csv(outfile, sep='\t', header=True, index=True, na_rep="NA")

def output_pheno_distributions(pcs, outgraphic, num_pcs):
    print("Outputting phenotype distributions to",outgraphic)
    # Determine plotting dimensions
    nplots = len(pcs.columns)
    ncol=num_pcs
    nrow = math.ceil(nplots / ncol)
    if (nrow-1) * ncol > nplots: nrow-=1    # Remove a row if unneeded

    # Set up figure
    fig = plt.figure(figsize=(ncol*4, nrow*3))
    grid = gridspec.GridSpec(nrows=nrow, ncols=ncol, hspace=0.5, wspace=0.5)

    # Plot each set of PCs as a histogram
    i=0
    for row in range(nrow):
        for col in range(ncol):
            if i >= nplots: break   # Since will have empty canvases at the end, break out once hit the end
            mydata = np.array(pcs.iloc[:,i])
            mydata = mydata[np.isfinite(mydata)]
            mytitle = pcs.columns[i]
            ax = fig.add_subplot(grid[row, col], xlabel="Distribution of values", ylabel="Count")
            ax.hist(mydata, bins=25, color=get_color(row))
            ax.set_title(mytitle, fontsize="xx-small", weight="bold")
            i+=1
    fig.savefig(outgraphic, dpi=100)


def output_scree_plots(eigens, infiles, outfile, num_eigens):
    print("Outputting scree plots to",outfile)
    # Determine plotting dimensions
    nplots = len(infiles)

    # Set up figure
    fig = plt.figure(figsize=(6, nplots*4))
    grid = gridspec.GridSpec(nrows=nplots, ncols=1, hspace=0.5, wspace=2)

    # titles
    left_trim = commonprefix(infiles)
    newnames = [re.sub(pattern="^" + left_trim, repl="", string=s) for s in infiles]

    # Plot each set of PCs as a histogram
    for row in range(nplots):
        ax = fig.add_subplot(grid[row, :], xlabel="PC #", ylabel="Fraction variance explained")
        ax.scatter(x=np.arange(num_eigens)+1, y=eigens[row][:num_eigens], c=get_color(row), s=40)
        ax.set_title(newnames[row])
    fig.savefig(outfile, dpi=100)

def get_color(i):
    colors=['b','g','r','c','m','y','k']
    value = i % 7
    return colors[value]

if __name__ == '__main__': main()