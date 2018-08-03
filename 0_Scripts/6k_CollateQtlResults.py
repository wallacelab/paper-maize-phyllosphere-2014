__author__ = 'jgwall'

import argparse
import numpy as np
import pandas as pd
import re

debug = False

def main():
    args = parse_args()
    print("Loading GWAS results from", len(args.infiles), "input files")
    results = [pd.read_table(i, index_col=0) for i in args.infiles]
    results = pd.concat(results)
    print("\t",len(results),"identified QTL loaded")
    results.to_csv(args.outprefix + ".all_found_qtn.txt", sep='\t', index=False)

    # Load simulated QTL
    qtl = pd.read_table(args.qtl)
    qtl['name'] = [q.replace(".","_") for q in qtl['name']]

    collated = collate_results(results, qtl, args.pval_cutoffs)
    collated.to_csv(args.outprefix + ".collated.txt", sep='\t', index=False)

    # Parse down to a simpler table with means and standard deviations
    parsed = collated[['name', 'num_qtl', 'h2', 'p_0.01']].copy()
    parsed['name'] = [re.sub(string=n, pattern="_perm.+", repl="") for n in parsed['name']]
    means = parsed.groupby(['h2', 'num_qtl']).mean()
    sd = parsed.groupby(['h2', 'num_qtl']).std()
    means.columns, sd.columns = ["mean"], ["sd"]
    pretty = pd.concat([means, sd], axis=1)
    pretty.to_csv(args.outprefix + ".pretty.txt", sep='\t')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*", help="List of files with parsed QTL GWAS results; outputs from 6j_ParseGwasResultsWithSimulatedQtl.py")
    parser.add_argument("-o", "--outprefix", help="Output file prefix")
    parser.add_argument("-q", "--qtl", help="File of simulated QTL")
    parser.add_argument("-p", "--pval-cutoffs", type=float, default=0.01, nargs="*", help="P-value cutoffs to calculate results at")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def collate_results(results, qtl, pval_cutoffs):
    print("Collating results at pvalue cutoffs of",pval_cutoffs)

    # Create dataframe to hold tallies
    tally=pd.DataFrame(index = qtl['name'], columns=pval_cutoffs)
    found_qtn = pd.DataFrame(index = qtl['name'], columns=pval_cutoffs)
    for trait in tally.index:
        for p in tally.columns:
            hits = (results['Trait'] == trait) & (results['empirical_pval'] <= p)
            tally.loc[trait, p] = np.sum(hits)
            found_qtn.loc[trait, p] = ",".join(list(results['Marker'].loc[hits]))

    # Sanity check to make sure things will be joined correctly
    name_mismatch = [a!=b for a,b in zip(qtl['name'], tally.index)]
    if(np.any(name_mismatch)): print("### WARNING! Names are mismatched and will not join properly!! ###")

    # Concatenate and return; also print a summary of results
    for p in tally.columns:
        print("\tFound", np.sum(tally.loc[:,p]),"QTN at p-value cutoff",p)
        qtl["p_" + str(p)] = list(tally.loc[:,p])
    for p in found_qtn.columns:
        qtl["found_qtn_" + str(p)] = list(found_qtn.loc[:,p])

    return(qtl)

if __name__ == '__main__': main()