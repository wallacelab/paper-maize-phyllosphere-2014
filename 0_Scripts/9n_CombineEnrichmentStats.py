__author__ = 'jgwall'

import argparse
import pandas as pd

debug = False


def main():
    args = parse_args()
    print("Combining enrichment stats for",len(args.aea),"AEA enrichment files and",len(args.fisher),"Fisher exact enrichment files")

    # Load data
    aea = pd.concat([load_aea(infile) for infile in args.aea])
    fisher = pd.concat([load_fisher(infile) for infile in args.fisher])
    annotations = sorted(set(aea.index) | set(fisher.index))  # Get unique labels and sort
    print("Identified", len(annotations), "unique annotations to search for")
    # aea = aea.set_index(['Annotation'])

    # Make output dataframe
    group, p_aea, p_fisher, p_bonf, p_fdr = list(), list(), list(), list(), list()
    for a in annotations:
        group.append(get_annotation(a, aea, "set"))
        p_aea.append(get_annotation(a, aea, "pval"))
        p_fisher.append(get_annotation(a, fisher, "pval"))
        p_bonf.append(get_annotation(a, fisher, "bonferroni"))
        p_fdr.append(get_annotation(a, fisher, "fdr"))
    combined = pd.DataFrame({"Annotation":annotations, 'Set':group, 'AEA p-value':p_aea, 'Fisher Exact p-value':p_fisher,
                             "Fisher Exact Bonferroni":p_bonf, 'Fisher Exact FDR':p_fdr})
    combined = combined[['Set','Fisher Exact p-value', 'Fisher Exact Bonferroni', 'Fisher Exact FDR', 'AEA p-value','Annotation']]

    # Filter out AEA terms with p-value of 1 and no Fisher Exact. (Outputting these seems like a bug in the AEA program)
    if args.remove_aea_pval_1:
        toremove = (combined['AEA p-value'] == 1) & ((combined['Fisher Exact p-value'] == "NA"))
        combined = combined.loc[~toremove, :]
        print("\tAfter removing AEA output with p-values of 1 (and no Fisher Exact score), have",len(combined),"annotations left")

    # Write out
    print("Writing combined file to", args.outfile)
    combined.to_csv(args.outfile, sep='\t', index=False)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--aea", nargs="*", help="List of parsed Annotation Enrichment Analysis output files")
    parser.add_argument("-f", "--fisher", nargs="*", help="List of parsed Fisher Exact test output files")
    parser.add_argument("--remove-aea-pval-1", default=False, action="store_true",
                        help="Whether to remove AEA enrichment files that have p-values of 1 and no corresponding Fisher Exact test")
    parser.add_argument("-o", "--outfile", help="Output file")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

# Load AEA results and cut down to just the needed info; also format the KO/COG group to be prettier
def load_aea(infile):
    data=pd.read_table(infile)
    formatted = pd.DataFrame({"Annotation":data['branch_name']})
    formatted['set'] = [parse_aea_signature(s) for s in data['sig_name']]
    formatted['pval'] = data['p_value']
    formatted = formatted.set_index('Annotation')
    return formatted

def parse_aea_signature(s):
    if s.startswith("ko_enrich"): return "KEGG Orthology"
    if s.startswith("cog_enrich"): return "COG"
    return s

def load_fisher(infile):
    data=pd.read_table(infile)
    data = data.set_index('term')
    return(data)

# Helper function to pull out a term's value from a different column in a data frame or NA if not there
def get_annotation(name, source, column):
    if name not in source.index:
        return "NA"
    return source.loc[name, column]

if __name__ == '__main__': main()