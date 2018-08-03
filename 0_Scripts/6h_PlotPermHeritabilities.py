__author__ = 'jgwall'

import argparse
import numpy as np
import matplotlib.pyplot as plt

debug = False


def main():
    args = parse_args()
    print("Plotting heritabilities between original and permutations")

    # Load data
    orig = get_herits(args.orig)
    perms = get_herits(args.perms)
    if len(orig) !=1:
        print("WARNING!! More than one heritability found in original file!")

    # Plot data
    plot_herits(orig, perms, args.outfile)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--orig", help="Original file with heritabilities for this trait")
    parser.add_argument("-p", "--perms", help="File with permuted heritability analyssi (from TASSEL)")
    parser.add_argument("-o", "--outfile", help="Output file to graph to")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def get_herits(infile):
    print("\tLoading heritabilities from",infile)
    IN = open(infile, "r")
    herits=list()

    # Process header
    header = IN.readline().strip().split('\t')
    markerID = header.index("Marker")
    geneticID = header.index("Genetic Var")
    residID = header.index("Residual Var")

    # Go through file 1 line at a time and capture heritabilities
    for line in IN:
        data = line.strip().split('\t')
        if data[markerID] != "None": continue   # Marker is "None" only for basic population/kinship term
        genetic_var, residual_var = float(data[geneticID]), float(data[residID])
        h2 = genetic_var / (genetic_var + residual_var)
        herits.append(h2)
    IN.close()
    print("\t\tExtracted",len(herits),"heritabilies")

    return(herits)


def plot_herits(orig, perms, outfile):
    print("Plotting heritabilities")
    np.random.seed(1)

    fig = plt.figure()
    ax = fig.add_subplot(111, title="Comparison of heritabilities", ylabel="h2")

    ax.scatter(x=0, y=orig, color="red")    # original
    offset = np.random.uniform(-0.25, 0.25, size=len(perms))
    ax.scatter(x=1 + offset, y=perms, color="blue")  # perms

    #Perms as scatters
    ax.violinplot(perms, positions=[2])  # Perms as violinplot

    # Prettify
    ax.set_xticks([0,1,2])
    ax.set_xticklabels(["Orig", "Perms", "Perms (violin)"])

    # Save
    fig.savefig(outfile)


if __name__ == '__main__': main()