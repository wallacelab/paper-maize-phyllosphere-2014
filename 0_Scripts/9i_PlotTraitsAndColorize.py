__author__ = 'jgwall'

import argparse
import matplotlib.pyplot as plt
import pandas as pd
import re
import sys


debug = False


def main():
    args = parse_args()

    # Load data
    print("Plotting traits against each other based on data in",args.infiles)
    data = [pd.read_table(infile, index_col=0).transpose() for infile in args.infiles]
    merged = pd.concat(data, axis=1, join='outer')

    # Confirm requested values are all in data
    for trait in [args.xvals, args.yvals, args.color]:
        if trait not in merged.columns:
            print("\tError! Could not locate trait",trait,"in supplied data!")
            sys.exit(1)
    print("\tAll requested traits found in data:",args.xvals, args.yvals, args.color)

    # Extract data
    xvals = merged.loc[:,args.xvals]
    yvals = merged.loc[:,args.yvals]
    colorkey = merged.loc[:,args.color]

    # Make scatterplot
    fig = plt.figure(figsize=(3,2.5))
    ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
    scatter = ax.scatter(x=xvals, y=yvals, c=colorkey, linewidths=0.5, s=40, cmap=args.cmap)
    colorbar = plt.colorbar(scatter, ax=ax, fraction=0.1, aspect=10, ticks=[])

    # Prettify
    ax.tick_params(labelbottom=False, labelleft=False, left=False, right=False, top=False, bottom=False)
    ax.set_xlabel(prettify(args.xvals), fontsize="small", weight='bold')
    ax.set_ylabel(prettify(args.yvals), fontsize="small", weight='bold')
    colorbar.set_label(prettify(args.color), rotation=270, weight='bold', fontsize='x-small', va='bottom')

    # Save
    fig.savefig(args.outfile, dpi=150)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*", help="Trait files to unify (traits in rows and samples in columns)")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-x", "--xvals", help="Trait for the X axis")
    parser.add_argument("-y", "--yvals", help="Trait for the Y axis")
    parser.add_argument("-c", "--color", help="Trait to use to colorize the scatterplot")
    parser.add_argument("--cmap", default='jet', help="Whcih color map to use for the color scale")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def prettify(label):
    label = re.sub(string=label, pattern=".+\\.", repl="")  # Strip higher-order taxonomy from taxa
    label = re.sub(string=label, pattern="_", repl=" ") # Change underscores back to spaces
    label=label.title() # Capitalization
    label = re.sub(string=label, pattern="Pc([0-9])$", repl="PC\\1") # Fix what title case does to "PC" in principal components
    return label



if __name__ == '__main__': main()