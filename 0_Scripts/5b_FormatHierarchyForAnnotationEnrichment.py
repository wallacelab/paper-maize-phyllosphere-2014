__author__ = 'jgwall'

import argparse
import networkx as nx
import numpy as np
import re

debug = False

# This script formats the metagenome hierarchy graph from 5a_AssembleMetagenomeHierarchy.py for Annotation Enrichment
#   Analysis (Glass & Grivan 2014, Sci Rep 4:4191, DOI: 10.1038/srep04191). For my purposes the lowest-level "genes" are
#   the individual KEGG or COG terms, and the higher-level hierarchies are the equivalent of the GO terms

# This script takes a graph file in gexf format (which should have 3 levels for COGs and 4 for KEGG pathways) and outputs
#    the following 3 files:
#       *.index.txt - A 3-column index file
#           Column 1 = ID# for each COG/KO term
#           Column 2 = Name of each term
#           Column 3 = Total # of higher-level terms they match to
#       *.annotations.txt - A 2-column file matching higher-level terms to individual COG/KO terms
#           Column 1 = The name of the higher-level term
#           Column 2 = Comma-separated list of which individual COG/KO terms match to this term (ie, are descendants of it)
#       *.partitions.txt - A 3-column file that contains the information necessary for the Annotation Enrichment program to build the graph
#           Column 1 - The name of the higher-level term (like in *.annotations.txt)
#           Column 2 - A comma-separated list of other higher-level terms this term is a descendant of. The first number is an index for
#               the term itself, so "Reproduction   17,257,199" means that Reproduction has ID 17 and is descended from terms with IDs 257 and 199


def main():
    args = parse_args()
    print("Converting metagenome graph from",args.infile,"into annotation enrichment input files")
    graph = nx.read_gexf(args.infile)

    lowest, higher, higher_key = get_graph_levels(graph, args.set)

    # Relabel nodes so work with AEA
    graph = nx.relabel_nodes(graph, higher_key)

    # Output key
    output_key(higher_key, args.outprefix + ".term_key.txt")

    # Make index file
    make_index_file(lowest, graph, args.outprefix + ".index.txt")

    # Make annotation file
    make_annotation_file(lowest, higher, graph, args.outprefix + ".annotations.txt")

    # Make partition file
    make_partition_file(higher, graph, args.outprefix + ".partitions.txt")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outprefix")
    parser.add_argument("-s", "--set", default="XX", help="Name for this set; used to make the pseudo-GO style annotation key (required by AEA)")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

# Take a graph and get all the lowest-level nodes and all the higher-level ones
def get_graph_levels(graph, name):
    # Determine graph depth
    lowest, higher = set(), set()
    for node in graph.nodes():
        ancestors = nx.ancestors(graph, node)
        descendants = nx.descendants(graph, node)

        # Sanity check to make sure is not a
        if len(ancestors) == 0 and len(descendants) == 0:
            print("\tWarning! Node",node,"has neither ancestors nor descendants and will be skipped!")
            continue

        if len(ancestors) > 0 and len(descendants) == 0:
            lowest.add(node)
        else: higher.add(node)

    print("\tLocated",len(lowest),"lowest-level nodes and",len(higher),"higher-level nodes")

    # Convert to dictionaries so save their unique ID #s
    lowest = make_id_dict(lowest)
    higher = make_id_dict(higher)

    # Make a key for the higher-level ones and create a replacement key with the new levels
    key = dict()
    name=name[:2] # Cut Set down to first 2 characters
    new_higher = dict()
    for h in higher:
        newname = name + ":" + str(higher[h]).zfill(7)
        key[h] = newname
        new_higher[newname] = higher[h]

    return lowest, new_higher, key

# Convert a set of nodes into a dictionary with unique IDs
def make_id_dict(nodes):
    result = dict()
    for n, i in zip(sorted(nodes), range(len(nodes))):
        result[n]=i
    return result


def make_index_file(nodes, graph, outfile):
    print("Making index file for",len(nodes),"nodes and writing to",outfile)
    OUT = open(outfile, "w")
    for node in sorted(nodes.keys()):
        id = nodes[node]
        count = len(nx.ancestors(graph, node))
        output = "\t".join([str(id), node, str(count)]) + "\n"
        OUT.write(output)
    OUT.close()


def make_annotation_file(lowest, higher, graph, outfile):
    print("Making annotation file and writing to",outfile)
    OUT = open(outfile, "w")
    for node in sorted(higher.keys()):
        # Filter descendants for only those that are in the lowest level (= individual COG/KO terms)
        descendants = nx.descendants(graph, node)
        tokeep = set()
        for d in descendants:
            if d in lowest: tokeep.add(d)
        # print(node)
        # print("Descendants:",descendants)
        # print("To keep:", tokeep)

        # Make comma-separated list ofindividual terms
        tokeep = sorted(tokeep)
        ids = [lowest[k] for k in tokeep]
        ids = ",".join([str(i) for i in ids])

        # Write out
        output = sanitize_name(node) + "\t" + ids + ",\n" # Example files have terminal comma, so added one
        OUT.write(output)
    OUT.close()


def make_partition_file(nodes, graph, outfile):
    print("Making partition file and writing to",outfile)
    OUT = open(outfile, "w")
    for node in sorted(nodes.keys()):
        # Get node IDs
        myID = nodes[node]
        ancestors = nx.ancestors(graph, node)
        other_ids = [nodes[a] for a in ancestors]
        ids = [myID] + other_ids

        # Format and output
        output = sanitize_name(node) + "\t" + ",".join([str(i) for i in ids]) + ",\n"  # Example files have terminal comma, so added one
        OUT.write(output)

    OUT.close()

# Remove potentially problematic characters
def sanitize_name(name):
    name = re.sub(pattern="[\\[\\] ]", repl="_", string=name)
    return name

def output_key(key, outfile):
    print("Outputting term key to",outfile)
    OUT = open(outfile, "w")
    OUT.write("original\tnew\n")
    for term in sorted(key.keys()):
        OUT.write(term + "\t" + key[term] + "\n")
    OUT.close()

# def validate_truncation(nodes):
#     nodes = sorted(nodes.keys())
#     truncated = [n[:10] for n in nodes]
#     unique = np.unique(truncated)
#     print("Truncation test: truncating higher-level nodes results in",len(unique),"from",len(nodes),"input labels")
#     if len(unique) != len(nodes):
#         print("\tWARNING!! This means the data is not suitable!!!")

if __name__ == '__main__': main()