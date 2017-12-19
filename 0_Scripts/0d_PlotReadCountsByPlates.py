__author__ = 'jgw87'
"""
Put my read counts into well format for visualizing any plate effects
"""

import argparse
import gzip
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import re
import sys
import matplotlib

def main():
    args = parse_args()
    reads=load_reads(args.infile)
    platekey = load_plate_key(args.keyfile)
    plot_plates(platekey, reads, args.outfile)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-k", "--keyfile")
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()
    return args

def load_reads(infile):
    print("Loading read counts from",infile)
    data=pd.read_csv(infile, sep='\t')
    reads=dict()
    for sample, count in zip(data['sample'], data['reads']):
         reads[sample]=count
    print("\tLoaded",len(reads),"sample counts")
    return reads


def load_plate_key(keyfile):
    print("Loading sample key from",keyfile)
    platekey=dict()
    data=pd.read_csv(keyfile, sep='\t')
    for i in range(len(data)):
        sample=data['Full Label'].iloc[i]
        plate=str(data['DNA_Plate'].iloc[i])
        row=data['Row (extrapolated)'].iloc[i]
        col=data['Column (extrapolated)'].iloc[i]
        if plate not in platekey:
            platekey[plate] = dict()
        if sample in platekey[plate]:
            print("Error! Sample",sample,"already loaded into key for",plate,"! Sample IDs should be unique!")
            sys.exit(1)
        platekey[plate][sample]=dict(plate=plate, row=row, col=col)
    for plate in platekey:
        print("\tLoaded",len(platekey[plate]),"samples for plate",plate)
    return platekey


def plot_plates(platekey, reads, outfile):
    print("Plotting plate heatmaps to",outfile)
    nrow, ncol=len(platekey), 1
    fig=plt.figure(figsize=(5, 4 * nrow))
    grid = gridspec.GridSpec(nrows=nrow, ncols=ncol, hspace=0.5, wspace=0.25)

    i=0
    for plate in sorted(platekey.keys()):
        print(plate)
        colvals, rowvals, depth = list(), list(), list()
        blanks=list()
        for sample in platekey[plate]:
            # print(sample)
            if sample not in reads: # Skip samples I don't have read depth for
                continue
            mydata = platekey[plate][sample]
            colvals.append(mydata['col'])
            rowvals.append(ord(mydata['row']))    # Convert character to int value for ease of plotting
            depth.append(float(reads[sample]))

            if "BLANK" in sample.upper():
                blanks.append(dict(col=mydata['col'], row=ord(mydata['row'])))

        # Check that actually have data to plot
        if len(colvals)==0: continue

        # Make figure axis and draw plots
        ax = fig.add_subplot(grid[i,:])
        ax.scatter(x=colvals, y=rowvals, c=np.log(depth), s=300, cmap=matplotlib.cm.coolwarm)

        # Add labels with 1000s of reads
        stringvals = [str(int(d/1000)) + "k" for d in depth]
        for x, y, s in zip(colvals, rowvals, stringvals):
            ax.text(x=x, y=y, s=s, fontsize="xx-small", horizontalalignment="center", verticalalignment="center")

        # Add title and X and Y axis labels
        ax.set_title(plate + " (log-scaled color)")
        ax.set_ylabel("Row")
        xlabels=np.unique(rowvals)
        ax.set_yticks(xlabels)
        ax.set_yticklabels([chr(x) for x in xlabels])
        ax.set_xlabel("Col")
        ax.set_xticks(np.unique(colvals))
        ax.invert_yaxis()

        # Mark blanks with red squares
        for b in blanks:
            ax.add_patch(matplotlib.patches.Rectangle(xy=[b['col']-0.5, b['row']-0.5], width=1, height=1,
                                                      facecolor="None", edgecolor="black", linewidth=3))

        i+=1
        # break

    fig.savefig(outfile, dpi=150, bbox_inches="tight")

if __name__ == "__main__":
    main()