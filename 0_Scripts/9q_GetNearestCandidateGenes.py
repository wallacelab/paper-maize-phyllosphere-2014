__author__ = 'jgwall'

import argparse
import numpy as np
import pandas as pd
import re

debug = False
# matplotlib.rcParams.update({"font.size":'10'})

def main():
    args = parse_args()
    print("Getting nearest candidate gene to each hit")

    # Load data
    data = pd.read_table(args.infile)
    genes=pd.read_table(args.gff, header=None)
    genes.columns=['seqname','source','feature','start','end','score','strand','frame','description']
    key=pd.read_table(args.key, index_col=0)

    # Subset hits to just those with right empirical p-value
    hits = data.loc[data['empirical_pval'] <= args.pval_cutoff].copy()  # Copy() to avoid Pandas complaining later
    print("\t",len(hits),"hits have empirical p-values <=",args.pval_cutoff)

    # Go through hits and find nearest n candidate genes
    candidates = get_candidate_genes(hits, genes, args.num_genes)
    print("Writing results")
    candidates.to_csv(args.outprefix + ".full.txt", sep='\t', index=False)

    # Make into a prettier version to coordinate hits with genes
    pretty = prettify_table(candidates, key)
    pretty.to_csv(args.outprefix + ".pretty.txt", sep='\t', index=False)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Input files of collected GWAS hits (best single hit per cluster, from step 9n_)")
    parser.add_argument("-o", "--outprefix", help="Output file prefix")
    parser.add_argument("-g", "--gff", help="GFF file of genes")
    parser.add_argument("-k", "--key", help="Keyfile of metagenome terms (to make them easier to read)")
    parser.add_argument("-n", "--num-genes", type=int, default=1, help="Number of closest genes to keep")
    parser.add_argument("-p", "--pval-cutoff", default=0.1, type=float, help="Cutoff for empirical p-values to be a 'hit'")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def get_candidate_genes(hits, genes, num_genes):
    candidates = list()
    for i in range(len(hits)):
        if debug and i > 5: break
        mychrom, mypos = str(hits['Chr'].iloc[i]), hits['Pos'].iloc[i]
        targets = genes.loc[genes['seqname'] == mychrom].copy()
        # print("\t",len(targets),"possible genes for",mychrom,":",mypos)

        # Determine if SNP is before, after, or inside of a gene
        before, after = targets['start'] > mypos, targets['end'] < mypos
        inside = (targets['start'] <= mypos) & (targets['end'] >= mypos)
        for j in range(len(targets)):
            if before.iloc[j] + after.iloc[j] + inside.iloc[j] != 1: print(
                "\tWarning! Locating SNP relative to genes gave >1 answer!")

        # Get distance from hit
        targets['distance'] = np.nan
        targets.loc[before, 'distance'] = targets['start'] - mypos
        targets.loc[after,  'distance'] = targets['end'] - mypos
        targets.loc[inside, 'distance'] = 0
        if np.any(np.isnan(targets['distance'])):
            print("\t\tWARNING! Found genes with nan distance to SNP")
        targets['abs_distance'] = np.abs(targets['distance'])
        # print(targets.head())

        # Sort based on absolute distance and take top n
        targets = targets.sort(columns='abs_distance')
        targets = targets.iloc[:num_genes, :]

        # Add info about hit
        targets['Trait'] = hits['Trait'].iloc[i]
        targets['Chr'] = mychrom
        targets['Pos'] = mypos
        targets['Marker'] = hits['Marker'].iloc[i]
        targets['p'] = hits['p'].iloc[i]
        targets['empirical_pval'] = hits['empirical_pval'].iloc[i]
        targets['trait_hit_cluster'] = hits['cluster'].iloc[i] if 'cluster' in hits.columns else "NA"
        targets['best_in_cluster'] = hits['best'].iloc[i] if 'best' in hits.columns else "NA"

        # Trim down and record
        output = targets.copy()
        cols_tokeep=['Trait', 'Marker', 'Chr', 'Pos', 'p', 'empirical_pval', 'trait_hit_cluster', 'best_in_cluster',
                     'feature', 'seqname', 'start', 'end', 'distance', 'description']
        output = output[cols_tokeep]
        candidates.append(output)

    candidates = pd.concat(candidates)
    return(candidates)


def prettify_table(candidates, key):
    print("Prettifying table for candidate gene output")
    pretty = candidates.sort(columns=['Trait','Chr','Pos']) # Sort
    print(pretty.head())
    pretty['Trait'] = clean_metagenome(pretty['Trait']) # Clean up names
    pretty = pretty.rename(columns={'p':'raw_pval', 'distance':'nearest_gene_distance'})


    # Get annotations for inferred metabolic functions
    annotations = key.loc[pretty["Trait"],:].copy()
    pretty['trait_description'] = np.array(annotations['description'])

    # Get annotations for candidate genes
    pretty['nearest_gene'] = get_geneID(pretty['description'])
    pretty['nearest_gene_description'] = get_gene_function(pretty['description'])

    # Drop extra columns and reorder
    pretty = pretty.drop(['feature', 'seqname', 'start', 'end','description'], axis='columns')
    pretty = pretty[['Trait','Marker','Chr','Pos','raw_pval','empirical_pval','trait_hit_cluster','best_in_cluster',
                     'trait_description','nearest_gene','nearest_gene_distance', 'nearest_gene_description']]
    return pretty

def clean_metagenome(traits):
    clean=list()
    for t in traits:
        t = re.sub(string=t, pattern="_[0-9]+$", repl="")
        t = re.sub(string=t, pattern="^log_", repl="")
        t = re.sub(string=t, pattern="_", repl=" ")
        t = re.sub(string=t, pattern="Meiosis   yeast", repl="Meiosis - yeast")
        t = re.sub(string=t, pattern="D Arginine and D ornithine", repl="D-Arginine and D-ornithine")
        clean.append(t)
    return(clean)


def get_geneID(descr):
    matches = [re.search(string=d, pattern="ID=gene:([^;]+);") for d in descr]
    names = [m.group(1) for m in matches]
    return names

def get_gene_function(descr):
    matches = [re.search(string=d, pattern="description=([^[]+) \[") for d in descr]
    names = [m.group(1) if m else "[no annotation]" for m in matches]
    names = [re.sub(string=n, pattern="%3B", repl=";") for n in names]
    names = [re.sub(string=n, pattern="%2C", repl=",") for n in names]
    return names

if __name__ == '__main__': main()