__author__ = 'jgwall'

import argparse
import numpy as np
import gzip

debug = False
metadata_cutoff=11

def main():
    args = parse_args()

    print("Randomly permuting genotypes in",args.infile,"with respect to taxa",args.n_perms,"times")
    np.random.seed(args.seed)
    IN = get_filehandle(args.infile, "rt")
    header=IN.readline().strip().split('\t')
    metadata=header[:metadata_cutoff]
    taxa=header[metadata_cutoff:]
    
    # Make permuted headers
    new_headers=list()
    OUTS=list()
    for i in range(args.n_perms):
        if debug and i > 5: break
        taxa = np.random.choice(taxa, size=len(taxa), replace=False)
        new_header = "\t".join(metadata + list(taxa)) + '\n'
        new_headers.append(new_header)

        outfile=args.outprefix + str(i) + args.outpostfix
        OUTS.append(get_filehandle(outfile, "wt"))

    # Write out the rest of the file
    genodata=IN.readlines()
    for myheader, OUT in zip(new_headers, OUTS):
        OUT.write(myheader)
        OUT.writelines(genodata)
        OUT.close()

    print("\tWrote permuted data to",len(OUTS),"output files")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outprefix", help="Part of output file to put before the permutation number")
    parser.add_argument("-p", "--outpostfix", default=".hmp.txt.gz", help="Part of output file to put after the permutation number")
    parser.add_argument("-n", "--n-perms", type=int, default=10, help="Number of permutations to make")
    parser.add_argument("-s", "--seed", type=int, default=1)
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def get_filehandle(file, mode):
    if file.endswith(".gz"):
        return gzip.open(file, mode)
    else: return open(file, mode)


if __name__ == '__main__': main()