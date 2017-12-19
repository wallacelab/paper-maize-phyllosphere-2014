__author__ = 'jgw87'
"""
Match the genotype names in a hapmap to sample names in my QIIME mapping file
"""

import argparse
import gzip
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import re

taxastart=11    # column in hapmap file where the taxa names start

def main():
    args = parse_args()
    map=pd.read_csv(args.mapfile, sep='\t')
    aliases = load_aliases(args.aliasfile)
    match = match_names(map, args.genofile, aliases)
    match.to_csv(args.outfile, sep='\t', index=False)
    if args.update_genofile:
        update_genofile(args.genofile, match, args.update_genofile)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mapfile")
    parser.add_argument("-g", "--genofile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-a", "--aliasfile", help="Two-column list of aliases to ehlp identify samples")
    parser.add_argument("-n", "--fix-names", help="Fix names to remove illegal characters")
    parser.add_argument("-u", "--update-genofile", required=False, help="(Optional) output file to restructure the hapmap")
    args = parser.parse_args()
    return args

def get_filehandle(file, mode, compression=5):  # Default to middle compression to go faster
    if file.endswith(".gz"):
        return gzip.open(file, mode, compression)
    else: return open(file, mode)

def match_names(map, genofile, aliases):
    print("Matching sample names to genotype name")
    # Matching genotype names
    match = map[['#SampleID',"Description"]].copy()
    match.columns = ["sample","name"]
    match["genotype"] = ""

    # Get header of genotype table
    IN = get_filehandle(genofile, "rt")
    header = IN.readline().strip().split('\t')
    header = header[taxastart:]
    taxa = [t.split(":")[0].upper() for t in header]   # Grab just the first part of the name & uppercase it
    #print(taxa)
    IN.close()

    # Go through and create a 3-column map of "this sample" = "maize line" = "this genotype"
    notfound = set()
    for i in range(len(match)):
        target = str(match['name'].iloc[i]).upper().replace(" ","")
        target = re.sub(string=target, pattern="GOODMAN-BUCKLER", repl="")
        # print(target,"yield alias",get_alias(target, aliases, taxa))
        if target == "NAN": continue    # Skip blanks
        if target in taxa:
            match['genotype'].iloc[i] = header[taxa.index(target)]
        elif get_alias(target, aliases, taxa) in taxa:
            match['genotype'].iloc[i] = header[taxa.index(get_alias(target, aliases, taxa))]
        else: notfound.add(target)

    print("\tUnable to find matches for", np.sum(match["genotype"]==""),"out of",len(match),"samples; Unique names are:", sorted(notfound))
    # print(match["name"].loc[match["genotype"]==""])

    return match

def update_genofile(infile, match, outfile):
    print("Updating hapmap to have sample names")
    print("\tLoading old genotypes")
    oldgenos = pd.read_csv(get_filehandle(infile, "rt"), sep='\t', low_memory=False)

    # Create new genotypes
    newgenos = oldgenos.iloc[:,:taxastart] # Copy over metadata; leave as a slice(?) to save memory
    for sample, genoname in zip(match['name'], match['genotype']):
        if (genoname in oldgenos.columns) and (sample not in newgenos.columns):
            newgenos[sample] = oldgenos[genoname]
    old_count = len(oldgenos.columns) - taxastart + 1
    new_count = len(newgenos.columns) - taxastart + 1
    print("\tNew genotype file has", new_count,"taxa versus the old",old_count,"taxa")

    # Write out
    print("\tWriting to",outfile)
    newgenos.to_csv(get_filehandle(outfile, "wt"), sep='\t', index=False)


def load_aliases(aliasfile):
    aliases = dict()
    if aliasfile is None: return aliases  # If no aliases, return an empty dictionary
    # Parse file
    print("Loading aliases from",aliasfile)
    IN = open(aliasfile, "r")
    header=IN.readline()
    for line in IN:
        v1, v2 = line.strip().split('\t')
        v1 = v1.upper()
        v2 = v2.upper()
        if v1 not in aliases:
            aliases[v1] = list()
        if v2 not in aliases:
            aliases[v2] = list()
        aliases[v1].append(v2)
        aliases[v2].append(v1)
    print("\tLoaded",len(aliases),"unique aliases (will be more than the number of lines in the file b/c stores in both directions)")
    # print(aliases)
    return aliases



# hack(ish) list of exceptions where simple rules don't match genotypes
exceptions = {


}
def get_alias(target, aliases, taxa):
    dummy_name = "TARGETNOTFOUND!!!!"
    target = target.upper()
    # If it isn't there, return a dummy sample name
    if target not in aliases:
        return dummy_name   # baiscally guaranteed to not be a name
    # If found and has only a single entry, return that
    if len(aliases[target]) == 1:
        return aliases[target][0]
    # If found and has multiple entries, return one that is in the taxa list (assuming it can be found)
    if len(aliases[target]) > 1:
        for t in aliases[target]:
            if t in taxa: return t
    # If everything up to here fails, return a dummy name
    return dummy_name

if __name__ == "__main__":
    main()