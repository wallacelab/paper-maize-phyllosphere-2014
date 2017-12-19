__author__ = 'jgwall'

import argparse
import pandas as pd

debug = False


def main():
    args = parse_args()

    # Load traits to keep
    traits = args.traitnames
    if args.traitfile:
        traits += [t.strip() for t in open(args.traitfile)]
    print("Filtering",args.infile,"for",len(traits),"traits")

    # Filter file
    data = pd.read_csv(args.infile, sep='\t', skiprows=1, header=None)
    oldtraits = data.iloc[1,:]
    traits = ["Taxon", "Taxa"] + traits # so keep the taxon column
    tokeep = [(o in traits) and (o not in args.exclude) for o in oldtraits] # Make true/false array of which ones to keep
    filtered = data.loc[:,tokeep]

    # Write out results
    print("\tWriting",len(filtered.columns),"traits to",args.outfile)
    open(args.outfile, "w").write("<Phenotype>\n")
    filtered.to_csv(args.outfile, sep='\t', index=False, mode="a", header=False, na_rep="NA")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-x", "--exclude", nargs="*", default=[], help='list of trait names to exclude')
    parser.add_argument("-t", "--traitfile", help="File of trait names to keep, 1 per line")
    parser.add_argument("-n", "--traitnames", nargs="*", default=[], help="List of trait names to keep")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


if __name__ == '__main__': main()