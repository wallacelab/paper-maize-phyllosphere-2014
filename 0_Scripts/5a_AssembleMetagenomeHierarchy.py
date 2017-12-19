__author__ = 'jgwall'

import argparse
import biom
import json
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

debug = False

# Notes: All the BIOM levels are length-1 numpy arrays, but structure is held in their format
#   Within each set of annotations, levels go from general (index 0) to specific (index 1 or 2 for COG or KO)

def main():
    args = parse_args()
    print("Making metagenome key for metagenome samples in",args.infile)

    table = biom.load_table(args.infile)
    names = table.ids(axis='observation')
    meta = table.metadata(axis='observation')
    if len(names) != len(meta) : print("\tWARNING! Length of observations and metadata do not match!!!")

    # Determine if is COGs or KOs
    keys = set()
    for m in meta:
        for k in m.keys():
            keys.add(k)
    group='unknown'
    if 'COG_Description' in keys: group='COG'
    elif 'KEGG_Description' in keys: group='KEGG'
    else: print("\tWARNING! Unknown set of metagenome annotations (neither COG nor KO)")

    # Make key components
    descriptions, categories, graph = parse_metadata(names, meta, group=group)

    # Assemble output table and write out
    print("\tAssembled key for",len(names),"IDs")
    output = pd.DataFrame({"id":names})
    output['description'] = descriptions
    output['unique_categories'] = [" | ".join(c) for c in categories]
    output.to_csv(args.outprefix + ".txt", sep='\t', index=False)

    # Write out graph of metagenome hierarchy
    print("\tMetagenome placed into a directed graph with",len(graph.nodes()),"nodes and",len(graph.edges()),"edges")
    nx.write_gexf(graph, args.outprefix + ".hierarchy.gexf")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outprefix")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def parse_metadata(names, meta, group):

    # Determine which names to use to access metadata
    if group == "COG":
        descr_name, cat_name = 'COG_Description', 'COG_Category'
    elif group == "KEGG":
        descr_name, cat_name = 'KEGG_Description', 'KEGG_Pathways'
    else: print("WARNING! unknown group",group)


    # Descriptions
    orig_descriptions = [m[descr_name] for m in meta]
    descriptions = list()
    for d in orig_descriptions:
        if len(d) != 1: print("WARNING! Description of length !=1 for",d)
        d = json.loads(d[0])

        # KO terms have multiple ones so have to join
        if type(d) == list:
            d = " | ".join(d)
        if type(d) != str: print("WARNING!  Non-string description for",d)
        descriptions.append(d)

    # Categories/Pathways, used to build up the hierarchy
    graph = nx.DiGraph()
    categories = list()
    for name, line in zip(names, meta):
        graph.add_node(name)

        # Parse out the annotations into a useful python data structure
        groups = line[cat_name]
        if len(groups) != 1: print("WARNING! Category of length != 1 for",groups)
        annotations = json.loads(groups[0])

        # Loop over the available annotation pathways (often just 1, but sometimes more) to build graph
        labels = set()
        for anno in annotations:
            # Sanity checks
            if group == "COG" and len(anno) != 2: print("WAWRNING! COG annotation",anno,"is not of expected length 2!")
            if group == "KEGG":
                if len(anno) == 1 and anno[0] == "None": pass   # Having "None" is okay
                elif len(anno) != 3: print("WARNING! KEGG annotation", anno, "is not of expected length 3!")

            # Rename annotations soas to preserve hierarchy in the face of duplicate names
            new_anno=anno.copy()
            for i in range(len(anno)):
                new_anno[i] = "->".join(anno[:i + 1])
            anno = new_anno

            for i in range(len(anno)):
                labels.add(anno[i]) # Add to set of labels for this COG/KO term
                graph.add_node(anno[i]) # Make a node in the ongoing graph



                # Add edge from previous level if this isn't the first level
                if i != 0:
                    graph.add_edge(anno[i-1], anno[i])

                # if this is the last level, add an edge from it to the individual COG/KO term
                if i == len(anno)-1 and anno[i] != "None":  # Don't want a bunch of things connecting to each other through "None"
                    graph.add_edge(anno[i], name)

        categories.append(sorted(labels))
    return descriptions, categories, graph


if __name__ == '__main__': main()