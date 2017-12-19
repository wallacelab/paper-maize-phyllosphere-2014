__author__ = 'jgwall'

import argparse
import matplotlib
import matplotlib.colors
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

debug = False
matplotlib.rcParams.update({'font.size': 8})

def main():
    args = parse_args()
    print("Plotting PCs")

    # Load data
    pcs=dict()
    for infile in args.infiles:
        name, data = load_pcs(infile)
        pcs[name]=data
    metrics = sorted(pcs.keys())
    key = pd.read_table(args.keyfile, index_col=0)

    # Plot out graphics
    fig = plt.figure(figsize=(15,6))
    grid = gridspec.GridSpec(nrows=len(args.categories), ncols=len(metrics)+1, hspace=.25, wspace=0.25)
    markerline=0.25
    markerscale=60
    panel_letter="A"
    for category in args.categories:
        row = args.categories.index(category)
        for metric in metrics:
            col = metrics.index(metric)
            mydata = pcs[metric]
            ax = fig.add_subplot(grid[row,col])
            pc1, pc2 = list(mydata.iloc[:,0]), list(mydata.iloc[:,1])
            colors, colorkey = get_colors(mydata.index, category, key)
            ax.scatter(pc1, pc2, c=colors, linewidths=markerline, s=markerscale)

            # Prettify axes
            title = (category + " (" + metric + ")").replace("_", " ").title()
            # ax.set_title(title, fontsize='large', fontweight='bold')
            ax.set_title(metric.title(), fontsize='large', fontweight='bold')
            ax.set_xlabel("PC1", fontweight="bold", size='medium')
            ax.set_ylabel("PC2", fontweight="bold", size='medium')
            ax.xaxis.set_ticks_position('none')  # bottom, top, both, or none
            ax.yaxis.set_ticks_position('none')  # left, right, both, or none
            ax.tick_params(labelleft='off', labelbottom='off')

            ax.text(-0.05, 1.09, panel_letter, transform=ax.transAxes, fontsize="x-large", fontweight='bold', va='top', ha='right')
            panel_letter = chr(ord(panel_letter) + 1)   # Increment panel letter to B, C, D, etc.

        # Make dummy plot as a legend for colors
        ax = fig.add_subplot(grid[row, len(metrics)])
        ax.axis('off')  # turn off border lines
        for label in get_order(colorkey.keys()):
            ax.scatter(0,0, c=colorkey[label], label=fancy_label(label), linewidths=markerline, s=markerscale)
        ax.set_xlim(1,2)    # So no points visible
        legend = ax.legend(loc='upper center', scatterpoints=1, frameon=False, title=fancy_label(category))
        leg_title = legend.get_title()
        leg_title.set_fontweight("bold")
        leg_title.set_fontsize("large")
        # ax.set_title(fancy_label(category), fontweight='bold')

        # Scoot axes over a little
        mypos = ax.get_position()
        mypos.x0 -= 0.07
        mypos.x1 -= 0.07
        ax.set_position(mypos)

    fig.savefig(args.outprefix + ".png", dpi=150, bbox_inches="tight")
    fig.savefig(args.outprefix + ".svg", dpi=150, bbox_inches="tight")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*", help="QIIME-created PCA files")
    parser.add_argument("-o", "--outprefix", help="Output file prefix")
    parser.add_argument("-k", "--keyfile", help="QIIME-formatted keyfile for sample data")
    parser.add_argument("-c", "--categories", nargs="*", help="Columns in keyfile to use for plotting")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


# Load and parse QIIME-formatted PC file
def load_pcs(infile):
    print("\tLoading PCS from",infile)
    IN = open(infile)
    while "Site\t" not in IN.readline(): pass   # Skip past header data

    # Read in
    samples, pcs = list(), list()
    for line in IN:
        if line == "\n": break  # Stop once hit bottom of PC block
        data=line.strip().split()
        samples.append(data[0])
        pcs.append(data[1:])

    # Convert to Pandas dataframe
    data = pd.DataFrame(pcs, index=samples)
    data.columns = ["PC" + str(i) for i in range(1, len(data.columns)+1)]

    # Get distance metric name
    name = re.sub(string=infile, pattern=".+/(.+)_pc.txt$", repl="\\1")
    return name, data

# Parse QIIME key to get colors out
def get_colors(samples, category, key):
    if category not in key.columns:
        print("\tWarning!",category,"not found in QIIME key")
        return "blue"

    # Get list of labels
    labels = np.array([key.loc[s, category] for s in samples])
    levels = sorted(set(labels))

    # Make color key
    np.random.seed(1)
    colornames = sorted(matplotlib.colors.cnames.keys())
    colorvals = np.random.choice(colornames, size=len(levels), replace=False)

    # Make colorkey
    colorkey = {level: value for level, value in zip(levels, colorvals)}

    # Manual override for specific ones I want to plot in specific ways
    if "KAK_7_8" in labels:
        colorkey={"KAK_7_8": "brown",
                  "KAK_9":   "peru",
                  "KAK_11":  "goldenrod",
                  "KAK_12": "darkseagreen",
                  "KAK_10": "lightseagreen",
                  "KAK_13":  "dodgerblue",
                  "KAK_14+": "mediumorchid"}
    elif "LMAN_8" in labels:
        colorkey = {"LMAN_8": "cornflowerblue",
                    "LMAN_26": "darkblue",
                    "LMAD_8": "gold",
                    "LMAD_26": "darkorange"}

    # Make list of colors
    colors = [colorkey[l] for l in labels]
    colors = matplotlib.colors.ColorConverter().to_rgba_array(colors, alpha=1)


    return colors, colorkey

# Custom order for legend
def get_order(labels):
    if "LMAD_8" in labels:
        return ["LMAD_8", "LMAN_8", "LMAD_26", "LMAN_26"]
    if "KAK_7_8" in labels:
        return ["KAK_7_8", "KAK_9","KAK_11","KAK_12","KAK_10","KAK_13","KAK_14+"]
    return sorted(labels)

# Fancy labels for legend
def fancy_label(label):

    # Day/night
    if "LMAD" in label: return re.sub(pattern="LMAD_(.+)", string=label, repl="Aug \\1 - Day")
    if "LMAN" in label: return re.sub(pattern="LMAN_(.+)", string=label, repl="Aug \\1 - Night")

    # Plates
    if label == "KAK_7_8": return "Batch 1 Plate 1"
    if label == "KAK_9": return "Batch 1 Plate 2"
    if label == "KAK_11": return "Batch 1 Plate 3"
    if label == "KAK_12": return "Batch 1 Plate 4"

    if label == "KAK_10": return "Batch 2 Plate 1"
    if label == "KAK_13": return "Batch 2 Plate 2"
    if label == "KAK_14+": return "Batch 2 Plate 3"

    # Category names
    if label=="dna_plate": return "RT-PCR Plate"
    if label == "tissue_date" : return "Sample time"

    # Default if nothing found
    return label

if __name__ == '__main__': main()