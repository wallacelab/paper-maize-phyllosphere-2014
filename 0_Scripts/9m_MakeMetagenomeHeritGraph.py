__author__ = 'jgwall'

import argparse
import networkx as nx
import matplotlib.pyplot as plt
import math
import pandas as pd
import re

debug = False

root_key={
    "Metabolism" : {'color': {'r': 0, 'g': 0, 'b': 255, 'a': 0}},
    "Unclassified" :    {'color': {'r': 255, 'g': 255, 'b': 0, 'a': 0}},
    "Cellular Processes"  : {'color': {'r': 0, 'g': 255, 'b': 0, 'a': 0}} ,
    "Environmental Information Processing" : {'color': {'r': 0, 'g': 255, 'b': 255, 'a': 0}},
    "Genetic Information Processing" : {'color': {'r': 255, 'g': 0, 'b': 255, 'a': 0}},
    "Human Diseases": {'color': {'r': 0, 'g': 0, 'b': 0, 'a': 0}}
}

def main():
    args = parse_args()
    print("Making subgraph of significant hits in",args.infile)

    # Load heritability and graph
    targets = load_herits(args.infile, args.max_herit, args.pval_cutoff)
    graph = nx.read_gexf(args.graphfile)

    # Find targets among renamed nodes with full hierarchy
    targets = find_targets(graph, targets)

    # Get sets of terminal significant ones and their ancestors
    ancestors = set()
    nodes = graph.nodes()
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

    # Prettify subgraph
    subgraph = prettify_subgraph(subgraph, targets)

    # Output
    #subgraph = subgraph.reverse()   #Flip so individual nodes have arrows toward their supercategory, not vice-versa
    nx.write_gexf(subgraph, args.outprefix + ".gexf", version="1.2draft")
    nx.draw_networkx(subgraph)
    plt.savefig(args.outprefix + ".png")

    # Output version without terminal leaves
    remove_terminal_leaves(subgraph)
    nx.write_gexf(subgraph, args.outprefix + ".trimmed.gexf", version="1.2draft")
    nx.draw_networkx(subgraph)
    plt.savefig(args.outprefix + "trimmed.png")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Parsed heritabilities from 5_AssembleMetagenomeHeritData.py")
    parser.add_argument("-g", "--graphfile", help="GEXF-formatted graph file of the annotation hierarchy")
    parser.add_argument("-f", "--fisher", help="P-value results for Fisher-exact test of term enrichment")
    parser.add_argument("-o", "--outprefix", help="Output file prefix")
    parser.add_argument("-p", "--pval-cutoff", type=float, default=1, help="P-value cutoff for including in the signature")
    parser.add_argument("-m", "--max-herit", type=float, default=1, help="Exclude heritability above this (b/c is probably a artifact)")
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


def find_targets(graph, targets):
    prettykey = {pretty_name(node):node for node in graph.nodes()}
    new_targets=set()
    for t in targets:
        if t in prettykey:
            new_targets.add(prettykey[t])
        else:
            print("Warning! Unable to find",t)
    return(new_targets)

# Make graph pretty
def prettify_subgraph(graph, hits):

    # Relabel nodes with shorter names
    label_key = {n:pretty_name(n) for n in graph.nodes()}   # Remove higher-level pathway info
    label_key = {n: remove_base_names(label_key[n]) for n in label_key} # Strip base COG/KO terms
    nx.set_node_attributes(graph, name="label", values=label_key)


    # Colorize nodes by their core region (metabolism, etc)
    for n in graph.nodes():
        found = 0
        colorkey = {'color': {'r': 125, 'g': 125, 'b': 125, 'a': 0}}    # Default gray
        # Colorize hits red
        if n in hits or pretty_name(n) in hits:
            colorkey = {'color': {'r': 255, 'g': 0, 'b': 0, 'a': 0}}    # Red if hit
            found+=1
        # Colorize other nodes according to their ancestral group
        else:
            ancestors = nx.ancestors(graph, n)
            for root in root_key:
                if n == root or root in ancestors:
                    colorkey = root_key[root]
                    found+=1
        if found != 1:
            print("!!WARNING!! Found",found,"root nodes for",n,"with ancestors",ancestors)
        graph.node[n]['viz'] = colorkey
    # print("\tFound",found,"hits that passed heritability")

    # Set size based on cumulative number of hits
    countkey = dict()
    for n in graph.nodes():
        desc = nx.descendants(graph, n) # Set of descendant nodes
        myhits = desc & hits
        mysize = (len(myhits) + 1) * 5
        graph.node[n]['viz']['size'] = math.sqrt(mysize)
        countkey[n] = len(myhits)
    nx.set_node_attributes(graph, name="hit_count", values=countkey)


    return graph

def pretty_name(node):
    if "->" in node:
        node = node.split("->")[-1]
    if node == "D-Arginine and D-ornithine metabolism":
        node = "D Arginine and D ornithine metabolism"
    return node

def remove_base_names(node):
    if node.startswith("K0") or node.startswith("K1") or node.startswith("COG"):
        node=""
    return node

def remove_terminal_leaves(graph):
    for node in graph.nodes():
        if (node.startswith("K0") or node.startswith("K1") or node.startswith("COG")) and len(nx.descendants(graph, node)) ==0:
            graph.remove_node(node)


if __name__ == '__main__': main()