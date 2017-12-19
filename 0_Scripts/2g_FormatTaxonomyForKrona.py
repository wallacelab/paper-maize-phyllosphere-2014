__author__ = 'jgw87'
"""
Take an MDS coordinate file and (text) otu abundance table from QIIME to color according to various OTUs
"""

import argparse
import pandas as pd
import re
import biom

def main():
    args = parse_args()
    print("Formatting taxonomy of OTUs in",args.infile,"and",args.biomfile,"for Krona")
    # Load taxonomy
    taxonomy = load_taxonomy(args.infile, args.include_otus)

    # Load OTU counts
    otus = biom.load_table(args.biomfile)
    #def identity(x, y): print(x,y)
    total_counts=otus.sum(axis="observation")
    otu_ids = otus.ids("observation")
    countdata = {o.lstrip("0"):c for o, c in zip(otu_ids, total_counts)}    # Have to left- strip zeros to get to work

    print("Writing results to",args.outfile)
    OUT = open(args.outfile, "w")
    for taxon in countdata:
        if taxon not in countdata:  print("\tWarning:",taxon,"not found in biom file")
        if taxon not in taxonomy:  print("\tWarning:", taxon, "not found in taxonomy file")
        OUT.write(str(countdata[taxon]) + "\t" + "\t".join(taxonomy[taxon]) + "\n")
    OUT.close()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-b", "--biomfile")
    parser.add_argument("-u", "--include-otus", default=False, action="store_true", help="Whether to include individual OTU labels")
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()
    return args

def load_taxonomy(infile, include_otus):
    data = pd.read_csv(infile, sep='\t', index_col=0, header=None, dtype=object)
    # print(data.head())
    data.columns = ["taxonomy","unknown1","unknown2"]
    taxonomy=dict()
    for i in range(len(data)):
        otu = str(data.index[i])
        taxonomy[otu] = taxonomy_from_string(data["taxonomy"].iloc[i], otu, include_otus)
        # print(otu,"has taxonomy",taxonomy[otu])
    return taxonomy

def taxonomy_from_string(string, otu, include_otus):
    notfound=""
    kingdom = re.search(string=string, pattern="k__([^;]*)").group(1) if "k__" in string else notfound
    phylum = re.search(string=string, pattern="p__([^;]*)").group(1) if "p__" in string else notfound
    class_ = re.search(string=string, pattern="c__([^;]*)").group(1) if "c__" in string else notfound
    order = re.search(string=string, pattern="o__([^;]*)").group(1) if "o__" in string else notfound
    family = re.search(string=string, pattern="f__([^;]*)").group(1) if "f__" in string else notfound
    genus = re.search(string=string, pattern="g__([^;]*)").group(1) if "g__" in string else notfound
    species = re.search(string=string, pattern="s__([^;]*)").group(1) if "s__" in string else notfound

    formatted=[kingdom]
    if phylum != notfound: formatted.append(phylum)
    if class_ != notfound: formatted.append(class_)
    if order != notfound:  formatted.append(order)
    if family != notfound: formatted.append(family)
    if genus != notfound:  formatted.append(genus)
    if species != notfound: formatted.append(species)

    if include_otus:
        formatted.append(otu)

    # print(string)
    # print(formatted,"\n")

    return formatted


if __name__ == "__main__":
    main()