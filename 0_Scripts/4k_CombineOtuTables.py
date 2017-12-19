__author__ = 'jgwall'

import argparse
import pandas as pd

debug = False


def main():
    args = parse_args()
    print("Combining",len(args.infiles),"input OTU tables into a single one")
    otus = [pd.read_csv(infile, sep='\t', index_col=0, skiprows=args.skip) for infile in args.infiles]
    combined=pd.concat(otus)
    
    if args.make_unique:
        print("Appending a unique digit to each OTU to make unique")
        new_index = [i + "_" + str(n) for i, n in zip(combined.index, range(len(combined.index)))]
        combined.index = new_index
    
    if args.min_nonzero:
        print("Removing OTUs with less than",args.min_nonzero,"nonzero counts")
        iszero = combined == 0
        iszero = iszero.sum(axis=1)
        nonzero_fract = 1 - (iszero/len(combined.columns))
        oldlength = len(combined)
        combined = combined.loc[nonzero_fract >= args.min_nonzero,:]
        print("\t",len(combined),"OTUs remain of an original",oldlength)
    
    combined.to_csv(args.outfile, sep='\t', index_label="#OTU")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-s", "--skip", type=int, default=0,  help="Number of lines at beginning of file to skip (if has extra header rows, for example)")
    parser.add_argument("-u", "--make-unique", default=False, action="store_true", help="Whether to make OTU names unique")
    parser.add_argument("-z", "--min-nonzero", type=float, help="Filter OTUs based on a minimal fraction of nonzero reads (optional)")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


if __name__ == '__main__': main()
