__author__ = 'jgw87'
"""
Take an MDS coordinate file and (text) otu abundance table from QIIME to color according to various OTUs
"""

import argparse
import pandas as pd
import re


def main():
    args = parse_args()
    print("Collapsing taxonomy of OTUs in",args.infile)
    taxonomy = load_taxonomy(args.taxonomyfile)
    otus = pd.read_csv(args.infile, sep='\t', skiprows=1, index_col=0)
    # print(otus)
    collapsed = collapse_taxa(otus, taxonomy)

    # Add prefix if specified
    if args.prefix:
        print("Adding",args.prefix,"as a prefix to all OTUs")
        collapsed.index = [args.prefix + "." + str(otu) for otu in collapsed.index]

    print("Writing results to",args.outfile)
    collapsed.to_csv(args.outfile, sep='\t', index_label="#OTU")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-t", "--taxonomyfile")
    parser.add_argument("-p", "--prefix", help="Prefix to add to all taxa names (including original OTUs)")
    args = parser.parse_args()
    return args

def load_taxonomy(infile):
    data = pd.read_csv(infile, sep='\t', index_col=0, header=None)
    # print(data.head())
    data.columns = ["taxonomy","unknown1","unknown2"]
    taxonomy=dict()
    for i in range(len(data)):
        otu = str(data.index[i])
        taxonomy[otu] = taxonomy_from_string(data["taxonomy"].iloc[i])
        # print(otu,"has taxonomy",taxonomy[otu])
    return taxonomy

def taxonomy_from_string(string):
    notfound=""
    kingdom = re.search(string=string, pattern="k__([^;]*)").group(1) if "k__" in string else notfound
    phylum = re.search(string=string, pattern="p__([^;]*)").group(1) if "p__" in string else notfound
    class_ = re.search(string=string, pattern="c__([^;]*)").group(1) if "c__" in string else notfound
    order = re.search(string=string, pattern="o__([^;]*)").group(1) if "o__" in string else notfound
    family = re.search(string=string, pattern="f__([^;]*)").group(1) if "f__" in string else notfound
    genus = re.search(string=string, pattern="g__([^;]*)").group(1) if "g__" in string else notfound
    species = re.search(string=string, pattern="s__([^;]*)").group(1) if "s__" in string else notfound

    mytaxa = dict()
    mytaxa["kingdom"] = kingdom
    mytaxa["phylum"] = ".".join([kingdom, phylum])
    mytaxa["class"] = ".".join([kingdom, phylum, class_])
    mytaxa["order"] = ".".join([kingdom, phylum, class_, order])
    mytaxa["family"] = ".".join([kingdom, phylum, class_, order, family])
    mytaxa["genus"] = ".".join([kingdom, phylum, class_, order, family, genus])
    mytaxa["species"] = ".".join([kingdom, phylum, class_, order, family, genus, species])
    return mytaxa

# Collapse OTU counts by taxa, summing everything within each division
def collapse_taxa(otus, taxonomy):

    # Make a copy to update OTU designations with phylogeny, too
    new_otus = otus.copy()
    new_otus.index = [taxonomy[str(name)]['species'] + "." + str(name) for name in new_otus.index]

    print("\tOriginal table has",len(otus),"taxa and",len(otus.columns),"samples")
    collapsed=list()
    collapsed.append(new_otus)  # add OTUs as original thing
    for level in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]:
        descriptions = [taxonomy[str(otu)][level] for otu in otus.index]
        new_frame = otus.groupby(descriptions).sum()
        print("\t",len(new_frame),"categories when collapsed by",level)
        collapsed.append(new_frame)
    collapsed = pd.concat(collapsed)
    print("Final table has",len(collapsed),"taxa and",len(collapsed.columns),"samples")
    return collapsed


if __name__ == "__main__":
    main()