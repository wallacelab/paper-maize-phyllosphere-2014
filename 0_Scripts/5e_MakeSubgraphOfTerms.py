__author__ = 'jgwall'

import argparse
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

debug = False


def main():
    args = parse_args()
    print("Making subgraph of significant hits in",args.infile)

    # Load heritability and graph
    targets = load_herits(args.infile, args.max_herit, args.pval_cutoff)
    graph = nx.read_gexf(args.graphfile)

    # Get sets of terminal significant ones and their ancestors
    ancestors = set()
    nodes = set(graph.nodes())
    added=0
    for t in targets:
        if t not in nodes: continue # Skip any hits that aren't targets
        added+=1
        for a in nx.ancestors(graph, t):
            ancestors.add(a)
    print("\tAdded ancestors of",added,"nodes to look for")

    # Make subgraph
    subgraph = graph.subgraph(targets | ancestors)
    print("Subgraph has",len(subgraph.nodes()),"nodes of the original",len(graph.nodes()))

    # Output
    subgraph = subgraph.reverse()   #Flip so individual nodes have arrows toward their supercategory, not vice-versa
    nx.write_gexf(subgraph, args.outprefix + ".gexf")
    nx.draw_networkx(subgraph)
    plt.savefig(args.outprefix + ".png")

    # Debugging stuff
    print("Longest path legnth is",nx.dag_longest_path_length(graph))
    print("Longest path is", nx.dag_longest_path(graph))

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Parsed heritabilities from 5_AssembleMetagenomeHeritData.py")
    parser.add_argument("-g", "--graphfile", help="GEXF-formatted graph file of the annotation hierarchy")
    parser.add_argument("-o", "--outprefix", help="Output file prefix")
    parser.add_argument("-p", "--pval-cutoff", type=float, default=1,
                        help="P-value cutoff for including in the signature")
    parser.add_argument("-m", "--max-herit", type=float, default=1,
                        help="Exclude heritability above this (b/c is probably a artifact)")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def load_herits(infile, max_herit, pval_cutoff):
    data = pd.read_table(infile, index_col=0)

    # Reformat trait names
    traits = data['trait']
    traits = [re.sub(pattern="^log_", string=t, repl="") for t in traits]
    traits = [re.sub(pattern="_[0-9]+$", string=t, repl="") for t in traits]
    traits = [re.sub(pattern="___", string=t, repl=" - ") for t in traits]
    traits = [re.sub(pattern="_", string=t, repl=" ") for t in traits]
    data['trait'] = traits

    # Identify the enriched, non-enriched, and bad traits based on heritability and empirical p-value
    good_pval = data['empirical_pval'] <= pval_cutoff
    good_herit = data['h2'] <= max_herit
    data['classification'] = 'unknown'
    data.loc[good_pval & good_herit, 'classification'] = "target"
    data.loc[~good_pval & good_herit, 'classification'] = "background"
    data.loc[~good_herit, 'classification'] = "error"

    # Subset
    targets = data['trait'].loc[data['classification'] == 'target']
    background = data['trait'].loc[data['classification'] == 'background']
    errors = data['trait'].loc[data['classification'] == 'error']

    # Status report and sanity check
    print("Of", len(data), "base traits, identified", len(targets), "as good/enriched", len(background),
          "as non-enriched/background, and",
          len(errors), "as bad due to too high heritability")
    if len(data) != (len(targets) + len(background) + len(errors)):
        print("\tWARNING!! These three categories should sum to the total but they don't!!!")

    # Convert to sets and return
    targetset = {t for t in targets}
    if len(targetset) != len(targets): print("WARNING! Duplicate target names detected!!")
    return targetset



if __name__ == '__main__': main()