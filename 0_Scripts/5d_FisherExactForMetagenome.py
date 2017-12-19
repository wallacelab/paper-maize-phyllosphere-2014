__author__ = 'jgwall'

import argparse
import networkx as nx
import numpy as np
import pandas as pd
import re
from scipy.stats import fisher_exact
from statsmodels.sandbox.stats.multicomp import multipletests

debug = False


def main():
    args = parse_args()
    print("Performing Fisher Exact test on", args.set,"data in",args.infile)

    # Load targets and background lowest-level annotations
    targets, background = load_herits(args.infile, args.set, args.max_herit, args.pval_cutoff, args.outprefix + ".trait_categories.txt")

    # Load graph and get dictionaries with sets of all ancestor terms for each term
    print("Identifying ancestor terms based on hierarchy in", args.graphfile)
    graph = nx.read_gexf(args.graphfile)
    targets = get_ancestors(targets, graph)
    background = get_ancestors(background, graph)

    # Get a list of all ancestor terms in the targets
    target_terms = set()
    for t in targets:
        for myterm in targets[t]:
            target_terms.add(myterm)
    target_terms = sorted(target_terms)
    print("\tIdentified",len(target_terms),"terms to search for potential enrichemnt")

    # Go through and get significance
    significance=list()
    for t in target_terms:
        pval = calc_fisher_exact(t, targets, background)
        significance.append(pval)

    # Summarize
    output = pd.DataFrame({'term':target_terms})
    output['pval'] = significance
    output = output.sort(['pval'])
    print(np.sum(output['pval'] < 0.01), "terms have uncorrected p-values below 0.01.")

    # Multiple testing correction
    pvals = np.array(output['pval'])
    bonferroni = multipletests(pvals, method='bonferroni')
    fdr = multipletests(pvals, method='fdr_bh')
    output['bonferroni'] = bonferroni[1]
    output['fdr'] = fdr[1]

    # Output
    outfile = args.outprefix + ".pvals.txt"
    print("Writing result to",outfile)
    output.to_csv(outfile, sep='\t', index=None)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Parsed heritabilities from 5_AssembleMetagenomeHeritData.py")
    parser.add_argument("-g", "--graphfile", help="GEXF-formatted graph file of the annotation hierarchy")
    parser.add_argument("-s", "--set", choices=['cog','ko'], help="Whether this is COG or KEGG terms ")
    parser.add_argument("-o", "--outprefix", help="Output file prefix")
    parser.add_argument("-p", "--pval-cutoff", type=float, default=1, help="P-value cutoff for including in the signature")
    parser.add_argument("-m", "--max-herit", type=float, default=1, help="Exclude heritability above this (b/c is probably a artifact)")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def load_herits(infile, set, max_herit, pval_cutoff, outfile):
    data = pd.read_table(infile, index_col=0)
    # Reformat trait names
    traits = data['trait']
    traits = [re.sub(pattern="^log_", string=t, repl="") for t in traits]
    traits = [re.sub(pattern="_[0-9]+$", string=t, repl="") for t in traits]
    data['trait'] = traits

    # Subset out based on set
    pattern = ""
    if set == "ko": pattern = "K"
    elif set == "cog": pattern = "COG"
    else: print("Warning! Unable to recognize set", set, "to determine which temrs to filter for")
    tokeep = [t.startswith(pattern) for t in data['trait']]
    data = data.loc[tokeep, :]
    print("\t", len(data), "traits identified in set", set)

    # Identify the enriched, non-enriched, and bad traits based on heritability and empirical p-value
    good_pval = data['empirical_pval'] <= pval_cutoff
    good_herit = data['h2'] <= max_herit
    data['classification'] = 'unknown'
    data.loc[good_pval & good_herit, 'classification'] = "target"
    data.loc[~good_pval & good_herit, 'classification'] = "background"
    data.loc[~good_herit, 'classification'] = "error"

    # Write out data table with classifications
    data.to_csv(outfile, sep='\t', index=None)

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
    backgroundset = {b for b in background}
    if len(targetset) != len(targets): print("WARNING! Duplicate target names detected!!")
    if len(backgroundset) != len(background): print("WARNING! Duplicate backgrond names detected!!")
    return targetset, backgroundset


def get_ancestors(terms, graph):
    termkey = dict()
    for t in terms:
        termkey[t] = {a for a in nx.ancestors(graph, t)}

    return termkey


def calc_fisher_exact(term, targets, background):
    target_hits = np.array([term in targets[t] for t in targets])
    background_hits = np.array([term in background[b] for b in background])
    total = len(targets) + len(background)

    # # Debug check
    # print(term)
    # for t, hit in zip(targets.keys(), target_hits):
    #     print(t, hit)

    # Get counts
    target_with = np.sum(target_hits)
    target_without = np.sum(~target_hits)
    background_with = np.sum(background_hits)
    background_without = np.sum(~background_hits)

    table = np.array([[target_with, target_without] , [background_with, background_without]])
    if np.sum(table) != total:
        print("WARNING! Fisher exact table for",term,"does not add total of",total,"terms!:\n",table)

    oddsratio, pval = fisher_exact(table, alternative="greater")
    return(pval)


if __name__ == '__main__': main()