__author__ = 'jgwall'

import argparse
import pandas as pd

debug = False


def main():
    args = parse_args()
    print("Updating format of Annotation Enrichment from",args.infile)
    data = pd.read_table(args.infile, header=None)
    data.columns = ["sig_number","branch_number","sig_name","branch_name","sig_annotations","branch_annotations",
                    "shared_annotations","p_value"]

    # Load key for branch names
    keydata = pd.read_table(args.keyfile)
    key=dict()
    for original, new in zip(keydata['original'], keydata['new']):
        key[new] = original

    # Replace names
    new_names = [key[n] for n in data['branch_name']]
    data['branch_name'] = new_names

    # Create new output, moving p-values to front
    output = data.drop(["sig_number","branch_number"], axis=1)
    pvals = output['p_value']
    output = output.drop("p_value", axis=1)
    output.insert(0, "p_value", pvals)
    output = output.sort("p_value")
    output.to_csv(args.outfile, sep='\t', index=False)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-k", "--keyfile", help="Keyfile connecting full COG/KEGG term to its shortened name for AEA ")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


if __name__ == '__main__': main()