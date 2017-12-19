__author__ = 'jgwall'

import argparse
import numpy as np
import pandas as pd
import re

debug = False


def main():
    args = parse_args()
    print("Combining", args.type, "heritability files from", len(args.small), "and", len(args.big),"input files")

    # Load data and concatenate into a single table
    small = [pd.read_table(s) for s in args.small]
    big = [pd.read_table(b) for b in args.big]
    small = pd.concat(small)
    big = pd.concat(big)
    print("Before joining, small perm set has",len(small),"entries and big perm set has",len(big),"entries")

    # Combine based on type of dataset
    if args.type == "otu":
        output = combine_otus(big, small)
    elif args.type == "metagenome":
        print("\tMetagenome key:",args.metagenome_key)
        print("\tTrait list:", args.traits)
        key = pd.read_table(args.metagenome_key) if args.metagenome_key else None
        traits = set()
        for traitfile in args.traits:
            for line in open(traitfile):
                traits.add(line.strip().replace(" ","_").replace("-","_"))
        print(len(traits),"unique traits loaded")
        output = combine_metagenomes(big, small, key, traits)
        # print(key)
    else:
        print("WARNING! No applicable methods for combining tables of type",args.type)

    # Remove NA heritabilities
    to_remove = np.isnan(output['h2'])
    print("\t",np.sum(to_remove),"traits to be removed because have NA heritabilities")
    output = output.loc[~to_remove,:]

    # Output combined tables
    print("Final output has",len(output),"rows of data")
    output.to_csv(args.outfile, sep='\t', index=False)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--small", nargs="*", help="Small herit file of initial permutations (gets lower priority)")
    parser.add_argument("-b", "--big", nargs="*", help="Big herit file of larger permutations (gets higher priority)")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-t", "--type", choices=['otu', 'metagenome'], help='What type of trait this is')
    parser.add_argument("-k", "--metagenome-key", help='3-column key of metagenome categories')
    parser.add_argument("--traits", nargs="*", help='Single-column file(s) of traits to filter for')
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def combine_otus(big, small):
    kingdom, phylum, class_, order, family, genus, species, otu = list(), list(), list(), list(), list(), list(), list(), list()
    log, h2, pval = list(), list(), list()

    loaded = set()
    for mytrait, myh2, mypval in zip(big['trait'], big['h2'], big['empirical_pval']):
        loaded.add(mytrait)
        log.append("Y" if mytrait.startswith("log_") else "")  # Indicate if was log-transformed
        h2.append(round(myh2, 3))
        pval.append(mypval)
        k, p, c, o, f, g, s, n = extract_taxa_data(mytrait)
        kingdom.append(k)
        phylum.append(p)
        class_.append(c)
        order.append(o)
        family.append(f)
        genus.append(g)
        species.append(s)
        otu.append(n)

    # Not elegant to do it twice, but faster to copy-paste than make a single function
    skipped=0

    for mytrait, myh2, mypval in zip(small['trait'], small['h2'], small['empirical_pval']):
        if mytrait in loaded:
            skipped+=1
            continue

        log.append("Y" if mytrait.startswith("log_") else "")  # Indicate if was log-transformed
        h2.append(round(myh2, 3))
        pval.append(mypval)
        k, p, c, o, f, g, s, n = extract_taxa_data(mytrait)
        kingdom.append(k)
        phylum.append(p)
        class_.append(c)
        order.append(o)
        family.append(f)
        genus.append(g)
        species.append(s)
        otu.append(n)

    print("\tSkipped",skipped,"traits that were in big perm file")
    output = pd.DataFrame(
        {'Log_Transformed': log, "Kingdom": kingdom, "Phylum": phylum, "Class": class_, "Order": order,
         "Family": family, "Genus": genus, "Species": species, "OTU":otu, "h2": h2, "Empirical P-value": pval})
    output = output[['Log_Transformed', "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU", "h2",
                     "Empirical P-value"]]
    return output


def extract_taxa_data(trait):
    trait = re.sub(string=trait, pattern="^log_", repl="")
    trait = re.sub(string=trait, pattern="^qiime\\.", repl="")
    trait = re.sub(string=trait, pattern="New\\.ReferenceOTU", repl="New_ReferenceOTU")

    taxonomy = trait.split(".")
    if len(taxonomy) > 8: print("Warning! More than 8 levels found for:\n\t", trait, "\n\t", taxonomy)

    # Add placeholders for taxa levels that stop before OTUs
    if len(taxonomy) < 8:
        taxonomy += ['-'] * (8 - len(taxonomy))

    # Go through and add placeholders for unknown levels
    for i in range(len(taxonomy)):
        if taxonomy[i] == "": taxonomy[i] = "unknown"

    # Fix New.ReferenceOTU so is back to QIIME standard
    if "New_ReferenceOTU" in taxonomy[7]:
        taxonomy[7] = re.sub(string=taxonomy[7], pattern="New_ReferenceOTU", repl="New.ReferenceOTU")

    return taxonomy


def combine_metagenomes(big, small, key, traits):
    loaded = set()
    bigdata, loaded = load_metagenome(big, key, traits, loaded)
    smalldata, loaded = load_metagenome(small, key, traits, loaded)
    bigdata['Number Permutations'] = 10000
    smalldata['Number Permutations'] = 100

    output = pd.concat([bigdata, smalldata])
    output = output.sort(columns=["Empirical P-value", "h2"], ascending=[True, False])
    return output

# Helper function to load metagenome data files, since can get quite complicated
def load_metagenome(data, key, traits, loaded):

    log, names, h2, pval = list(), list(), list(), list()
    descriptions, categories = list(), list()
    to_expand = set(key['id']) if key is not None else None
    desc_key = np.array(key['description'])
    cat_key = np.array(key['unique_categories'])

    # Loop over traits
    for mytrait, myh2, mypval in zip(data['trait'], data['h2'], data['empirical_pval']):
        if mytrait in loaded: continue
        loaded.add(mytrait)

        # Extract information from heritability file
        name = re.sub(pattern="_[0-9]+$", string=mytrait, repl="")  # Remove terminal tag I added to make unique
        name = re.sub(pattern="^log_", string=name, repl="")

        if name not in traits:
            # print(name,"skipped because not in supplied traits")
            continue

        # Link to metagenome key file
        description, category = "unknown", "unknown"
        if (to_expand is not None) and (name in to_expand):
            target = np.where(key["id"] == name)[0]
            if len(target) > 1: print("WARNING! More than one key entry for", name)
            target = target[0]
            # print(name,"at",target)
            description = desc_key[target].replace("|", ";")
            category = cat_key[target].replace("|", ";")
            # print(key['category'].iloc[target])

        # Add to list
        log.append("Y" if mytrait.startswith("log_") else "")  # Indicate if was log-transformed
        h2.append(round(myh2, 3))
        pval.append(mypval)
        names.append(name)
        descriptions.append(description)
        categories.append(category)

    result = pd.DataFrame(
        {'Log_Transformed': log, "Name": names, "h2": h2, "Empirical P-value": pval, "Description":descriptions, "Category":categories})
    result = result[['Log_Transformed', "Name", "h2", "Empirical P-value", "Description", "Category"]]
    return result, loaded

if __name__ == '__main__': main()
