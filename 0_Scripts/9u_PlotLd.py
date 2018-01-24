__author__ = 'jgwall'

import argparse
import collections
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import re

debug = False
gff_header = ['chrom','source','type','start','stop','score','strand','frame','annotations']

def main():
    args = parse_args()
    print("Plotting LD from",args.infile,"with windows of", args.big_winsize,"and",args.small_winsize,"bp.")
    if args.small_winsize > args.big_winsize:
        print("\tWARNING! Second window should be smaller than the first; plot may not appear correct")

    # Load data and simplify
    ld = pd.read_table(args.infile, low_memory=False)
    genes = pd.read_table(args.gff, low_memory=False, comment="#", header=None, names=gff_header)
    ld['Locus1'] = [str(l) for l in ld['Locus1']]
    ld['Locus2'] = [str(l) for l in ld['Locus2']]
    target_chrom, target_pos = find_targets(ld)

    # Filter down to just target sites
    big_ld = filter_ld(ld, target_chrom, target_pos, args.big_winsize)
    small_ld = filter_ld(big_ld, target_chrom, target_pos, args.small_winsize)
    big_genes = filter_genes(genes, target_chrom, target_pos, args.big_winsize)
    small_genes = filter_genes(big_genes, target_chrom, target_pos, args.small_winsize)

    # write out subsetted data (mostly for error checking)
    big_ld.to_csv(args.outprefix + ".big_ld.txt", sep='\t', index=None)
    small_ld.to_csv(args.outprefix + ".small_ld.txt", sep='\t', index=None)

    # Set up figure
    fig = plt.figure(figsize=(8,3.25))
    grid = gridspec.GridSpec(nrows=100, ncols=100)
    ax_big = fig.add_subplot(grid[2:55, :], title="") #LD at " + str(target_chrom) + ":" + str(target_pos))
    ax_small = fig.add_subplot(grid[65:90, :], title="")

    # Plot LD
    plot_ld(ax_big, big_ld, target_pos)
    plot_ld(ax_small, small_ld, target_pos)

    # Clean up axes
    clean_axes(ax_big, target_pos, args.big_winsize)
    clean_axes(ax_small, target_pos, args.small_winsize)

    # Plot genes & save to file
    mygenes = plot_genes(ax_small, small_genes)
    mygenes.to_csv(args.outprefix + ".genes.gff", sep='\t', index=None)

    # Axes-specific stuff
    ax_big.spines['top'].set_visible(False)
    ax_big.spines['right'].set_visible(False)
    ax_small.spines['top'].set_visible(False)
    for spine in ['left','right','top','bottom']:
        ax_small.spines[spine].set_linewidth(0.5)
    ax_small.set_xlabel("Position (Megabases)", weight='bold')  # X-axis label
    plt.figtext(x=0.05, y=0.5, s="Linkage Disequilibrium (r$^2$)", weight='bold', rotation=90, verticalalignment='center')  #Y-axis label
    ax_small.axvline(target_pos, color='darkblue', zorder=-1, linewidth=0.5, linestyle="dashed")    # Line on small axes showing target
    add_inset_box(ax_big, ax_small, target_pos, args.small_winsize)
    ax_small.tick_params(axis='x', pad=10)
    ax_big.text(x=0.02, y=0.94, s="Position " + str(target_chrom) + ":" + str(target_pos), transform=ax_big.transAxes, weight='bold', size='x-small')

    # Save figure
    fig.savefig(args.outprefix + ".png", dpi=150)
    fig.savefig(args.outprefix + ".svg", dpi=600)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outprefix")
    parser.add_argument("-g", "--gff", help="GFF-formatted file of genes, ideally including exons, too")
    parser.add_argument("-b", "--big-winsize", type=int, default=10000000, help="Genomic widnow for big plot (in bp)")
    parser.add_argument("-s", "--small-winsize", type=int, default=100000, help="Genomic window for small plot (in bp)")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

# Function to find target site and chrom because TASSEL changes which column they're in
# Takes advantage of fact that these files are site-by-all, so the target site is the most common by far
def find_targets(ld):
    chroms = list(ld['Locus1']) + list(ld['Locus2'])
    poses = list(ld['Position1']) + list(ld['Position2'])

    mychrom, chromcount = collections.Counter(chroms).most_common(1)[0]
    mypos, poscount = collections.Counter(poses).most_common(1)[0]

    # Sanity checking
    if chromcount < len(chroms)/2: print("\tWarning! target chromosome",mychrom,"present in less than half the lines of data (",chromcount,"of",len(chroms),")")
    if poscount < len(poses) / 2: print("\tWarning! target position",mypos,"present in less than half the lines of data (",poscount,"of",len(poses),")")

    print("\tTarget is",mychrom,":",mypos,"; counts are",chromcount,"and",poscount,"out of",len(ld),"lines of data, respectively")
    return mychrom, mypos


def filter_ld(ld, target_chrom, target_pos, winsize):
    right_chrom = (np.array(ld['Locus1']) == target_chrom) & (np.array(ld['Locus2']) == target_chrom)
    right_pos1 = np.abs(np.array(ld['Position1']) - target_pos) <= winsize
    right_pos2 = np.abs(np.array(ld['Position2']) - target_pos) <= winsize
    right_pos = right_pos1 & right_pos2
    myld = ld.loc[np.array(right_chrom) & np.array(right_pos), :]
    print("\tIdentified", len(myld), "sites out of", len(ld), "to be included in the plot")
    return myld



def filter_genes(genes, target_chrom, target_pos, winsize):
    right_chrom = np.array(genes['chrom']) == target_chrom
    right_start = np.abs(np.array(genes['start']) - target_pos) <= winsize
    right_stop = np.abs(np.array(genes['stop']) - target_pos) <= winsize
    right_pos = right_start | right_stop
    mygenes = genes.loc[np.array(right_chrom) & np.array(right_pos), :]
    print("\tIdentified", len(mygenes), "gene features out of", len(genes), "to be included in the plot")
    return mygenes



def plot_ld(ax, ld, target_pos):
    xvals = list()
    for i in range(len(ld)):
        # Keep the x-value that is NOT the target
        if ld['Position1'].iloc[i] == target_pos:
            xvals.append(ld['Position2'].iloc[i])
        else:
            xvals.append(ld['Position1'].iloc[i])
    myscatter = ax.scatter(xvals, ld['R^2'], linewidth=0, color='darkred', alpha=0.5)
    myscatter.set_rasterized(True)
    ax.scatter(target_pos, 1, linewidth=1, color='red', edgecolors='black', clip_on=False)  # Actual SNP being interrogated


def plot_genes(ax, genes):
    y = -0.1   # Where to draw lines
    tokeep = list()
    for i in range(len(genes)):
        myline =  None
        mytype, mystart, mystop = genes['type'].iloc[i], genes['start'].iloc[i], genes['stop'].iloc[i]
        if mytype == 'gene' or mytype=='exon':
            tokeep.append(i)
            myline = Line2D(xdata=[mystart, mystop], ydata=[y, y], color='darkgreen', linewidth=2.5, solid_capstyle='butt', clip_on=False)
            # gene_name = re.search(pattern="ID=gene:([^;]+)", string=genes['annotations'].iloc[i]).group(1)
            # ax.text(x = (mystart + mystop)/2, y=y+0.04, s=gene_name, color='darkgreen', size='xx-small', weight='bold', ha='center')
        if mytype == 'exon':
            myline.set_linewidth(6)
        if myline:
            ax.add_line(myline)

    keepers = genes.iloc[tokeep,:]
    return keepers
    # ax.set_ylim(bottom = y - 0.05, top=1.05)





# Cleanup common to both plots
def clean_axes(ax, target_pos, winsize):
    # Title & axis labels
    ax.set_title(ax.get_title(), weight='bold')

    # Ticks and spines (=border lines)
    ax.xaxis.set_ticks_position('bottom')  # bottom, top, both, or none
    ax.yaxis.set_ticks_position('left')  # left, right, both, or none
    ax.tick_params(labeltop='off', labelright='off')


    # Axis limits
    ax.set_ylim(bottom=0, top=1.05)
    ax.set_xlim(left = target_pos - winsize*1.1, right= target_pos + winsize*1.1)

    # Deal with tick labels
    # x-axis
    xticklabels = list()
    for tick in ax.get_xticks():
        tick = round(tick / 1e6, ndigits=2)
        xticklabels.append("{:.2f}".format(tick))
    ax.set_xticklabels(xticklabels, weight='bold', size='x-small')

    # y-axis
    yticks, yticklabels = list(), list()
    for tick in ax.get_yticks():
        if tick <= 1:
            yticklabels.append(str(tick))
            yticks.append(tick)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, weight='bold', size='x-small')




# Add an inset box connecting the two graphics
def add_inset_box(ax_big, ax_small, target_pos, small_winsize):
    # Box around region of interest in big plot
    big_ylim=ax_big.get_ylim()
    inset = patches.Rectangle((target_pos-small_winsize, big_ylim[0]),  width=small_winsize*2, height=big_ylim[1]-big_ylim[0],
                              color='blue', fill=False, linewidth=0.5)
    ax_big.add_patch(inset)

    # Get lines connecting inset to zoom axes
    small_xlim, small_ylim = ax_small.get_xlim(), ax_small.get_ylim()


    rightend = ax_small.transData.transform([small_xlim[1], small_ylim[1]])

    # Left connector line
    leftstart = [target_pos - small_winsize, big_ylim[0]]
    leftend = ax_small.transData.transform([small_xlim[0], small_ylim[1]])
    leftend = ax_big.transData.inverted().transform(leftend)    # convert back to data coordinates
    left_connect = Line2D(xdata=[leftstart[0], leftend[0]], ydata=[leftstart[1], leftend[1]], clip_on=False, linewidth=0.5)
    ax_big.add_line(left_connect)

    # Right connector line
    rightstart = [target_pos + small_winsize, big_ylim[0]]
    rightend = ax_small.transData.transform([small_xlim[1], small_ylim[1]])
    rightend = ax_big.transData.inverted().transform(rightend)    # convert back to data coordinates
    right_connect = Line2D(xdata=[rightstart[0], rightend[0]], ydata=[rightstart[1], rightend[1]], clip_on=False, linewidth=0.5)
    ax_big.add_line(right_connect)

    rightstart = [target_pos + small_winsize, big_ylim[0]]



if __name__ == '__main__': main()