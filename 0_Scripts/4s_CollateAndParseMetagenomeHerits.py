__author__ = 'jgwall'

import argparse
import pandas as pd
import re

debug = False


def main():
    args = parse_args()
    print("Parsing metagenome results from",len(args.infiles), "input files")

    data=load_data(args.infiles, args.keep_flowering_time)
    traitkey = load_trait_key(args.cog_names, args.cog_functions, args.kegg)
    output = add_functions(data, traitkey)
    print("Writing output to",args.outfile)
    output.to_csv(args.outfile, index=False, sep='\t')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*", help="List of input files from metagenome heritability analysis")
    parser.add_argument("-o", "--outfile", help="Combined output file of heritable metagenome and functions")
    parser.add_argument("--cog-names", help="Filepath to COG names file from COG database")
    parser.add_argument("--cog-functions", help="Filepath to COG functions file from COG database")
    parser.add_argument("--kegg", help="Filepath to KEGG orthology file")
    parser.add_argument("--keep-flowering-time", default=False, action="store_true", help="Keep flowering time in the dataset")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


# Load data and remove flowering time if required
def load_data(infiles, keep_flowering_time):
    data = [pd.read_csv(i, sep='\t') for i in infiles]
    data=pd.concat(data)
    if keep_flowering_time:
       print("\tLoaded", len(data.index), "traits")
    else:
       data = data.loc[data['trait'] != 'flowering_time',:]
       print("\tLoaded", len(data.index), "traits after removing flowering time")
    return data


# Load explanations for COG and KEGG functionality
def load_trait_key(cog_names, cog_functions, kegg_annotations):
    print("Loading COG names and functions")
    traitkey = dict()

    # Functions first
    funcs = pd.read_csv(cog_functions, sep='\t')
    functionkey = {code:fun for code,fun in zip(funcs['# Code'], funcs['Name'])}
    for func in functionkey.values():   # Load higher-level functions into appropriate format for matching
        traitkey[func] = dict()
        traitkey[func]['function'] = func
        traitkey[func]['name'] = func

    # Now actual cog names
    cogs = pd.read_csv(cog_names, sep='\t')
    for cog, func, name in zip(cogs['# COG'], cogs['func'], cogs['name']):
        if cog in traitkey:  print("\tWARNING! COG",cog,"loaded more than once!")
        traitkey[cog] = dict()
        my_functions = [functionkey[f] for f in func]   # General class of COG, like "lipid transport and metabolism, etc; can have >1
        traitkey[cog]['function'] = ";".join(my_functions)
        traitkey[cog]['name'] = name

    # Now KEGG annotations
    kegg = pd.read_csv(kegg_annotations, sep='\t')
    for ko, func, name in zip(kegg['ko'], kegg['func'], kegg['name']):
        traitkey[ko] = dict()
        traitkey[ko]['name'] = name
        traitkey[ko]['function'] = func

    print("\tLoaded",len(traitkey),"traits to match")
    return traitkey


# Add functions to the heritabiltiy file
def add_functions(data, traitkey):
    print("Adding functions to heritability for analysis")
    names, functions = list(), list()
    for trait in data['trait']:
        oldtrait=trait
        trait = re.sub(string=trait, pattern="^log_", repl="")  # Remove log_ prefix for any ones that were log-transformed
        trait = re.sub(string=trait, pattern="_[0-9]+$", repl="")  # Remove unique identifier I added at an earlier step
        # print(oldtrait,"becomes",trait)
        if trait in traitkey:
            names.append(traitkey[trait]['name'])
            functions.append(traitkey[trait]['function'])
        else:
            print("\tWarning: Unable to locate function",trait)
            names.append('unknown')
            functions.append('unknown')
    data['function'] = functions
    data['name'] = names
    return(data)


if __name__ == '__main__': main()