__author__ = 'jgw87'
"""
Take a TASSEL taxasummary file and find all the taxa that share a primary name, then choose the one that has the best
coverage and output it to a text file. Also output a log to report some of the merged ones are _not_ the highest coverage
"""

import argparse
import gzip
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import re


def main():
    args = parse_args()
    print("Choosing the best taxa from",args.taxafile)
    coveragekey = load_taxa_summary(args.taxafile)
    good_taxa = choose_best_taxa(coveragekey)
    
    print("\tWriting full report to",args.outfile)
    good_taxa.to_csv(args.outfile, sep='\t', index=False)
    
    if args.bestfile:
        print("\tWriting best taxa to",args.bestfile)
        best=list(good_taxa['best'])
        OUT = open(args.bestfile, 'w')
        [OUT.write(b + '\n') for b in best]
        OUT.close()
        

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--taxafile")
    parser.add_argument("-o", "--outfile", help="Output file of entire tabulation of which taxon is the best")
    parser.add_argument("-b", "--bestfile", help="Outfile file with just the best taxa name (for downstream filtering)")
    args = parser.parse_args()
    return args

def load_taxa_summary(infile):
    taxasummary = pd.read_csv(infile, sep='\t')
    names = np.array(taxasummary['Taxa Name'])
    coverage = 1-np.array(taxasummary['Proportion Missing'])
    # Return as dictionary of name->coverage
    coveragekey=dict()
    for n, c in zip(names, coverage):
        coveragekey[n] = c
    return coveragekey

def choose_best_taxa(coveragekey):

    # Make a dictionary linking the base name to each GBS name
    namegroups = dict()
    for taxon in coveragekey:
        name = taxon.split(":")[0]
        if name not in namegroups: namegroups[name]=list()
        namegroups[name].append(taxon)
    print("\tFound",len(namegroups),"unique taxa out of",len(coveragekey),"unique GBS identifiers")

    # Go through taxa sets and choose the best from each
    good_taxa = dict(group=list(), best=list(), merged=list())
    for group in namegroups:
        best_taxon, merged_taxon = None, None
        # If only one taxon in the group, keep that and move on
        if len(namegroups[group]) == 1 :
            best_taxon = merged_taxon = namegroups[group][0]
        else:
            best_taxon, merged_taxon = choose_taxon(namegroups[group], coveragekey)
        good_taxa['group'].append(group)
        good_taxa['best'].append(best_taxon)
        good_taxa['merged'].append(merged_taxon)
    good_taxa['best_coverage'] = [coveragekey[t] for t in good_taxa['best']]
    good_taxa['merged_coverage'] = [coveragekey[t] for t in good_taxa['merged']]

    # Format things into DataFrame
    for key in good_taxa:
        good_taxa[key] = pd.Series(good_taxa[key])
    good_taxa=pd.DataFrame(good_taxa)
    good_taxa =good_taxa[['group','merged','merged_coverage','best','best_coverage']]
    good_taxa=good_taxa.sort(columns="group")
    good_taxa["difference"] = np.array(good_taxa["best_coverage"]) - np.array(good_taxa["merged_coverage"])

    # Check how many "best" and "merged" don't match
    flipped, flipnames = 0, list()
    for i in range(len(good_taxa)):
        if (good_taxa['best'].iloc[i] != good_taxa['merged'].iloc[i]) and \
            (good_taxa['best_coverage'].iloc[i] > good_taxa['merged_coverage'].iloc[i]):
            flipped+=1
            flipnames.append(good_taxa['group'].iloc[i])
    print("Found",flipped,"instances of the merged sample having lower coverage than the best one:", flipnames)
    print("\tThe biggest difference is",np.max(good_taxa['difference']))

    return good_taxa




def choose_taxon(group, coveragekey):
    best = group[0]
    merged = list()
    for taxon in group:
        if coveragekey[taxon] > coveragekey[best]:
            best = taxon
        if ":MRG:" in taxon:
            merged.append(taxon)

    if len(merged) > 1:
        depth = [int(taxon.split(":")[2]) for taxon in merged]
        merged = np.array(merged)[np.argsort(depth)[::-1]]  # Sort by decreasing depth
        print("WARNING! More than one merged taxon found:", merged)
    if len(merged) == 0: merged.append(best) # If didn't find any merged taxa, default to the best one
    return best, merged[0]


# def logit(LOG, message):
#     if LOG:
#         LOG.write(message + "\n")

if __name__ == "__main__":
    main()