__author__ = 'jgwall'

import argparse

debug = False


def main():
    args = parse_args()
    key = load_key(args.keyfile)

    update_blups(args.infile, key, args.outfile)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-k", "--keyfile", help="Key file matching line names to genotypes; from 4p_MatchGenoNamesToSamples.py")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def load_key(keyfile):
    print("Loading genotype key from",keyfile)
    IN = open(keyfile)
    header=IN.readline().strip().split('\t')
    nameID, genoID = header.index('name'), header.index('genotype')

    key=dict()
    for line in IN:
        data=line.strip('\n').split('\t')   # Only strip newlines; if strip everything, lose tabs that separate blank fields
        myname, mygeno = data[nameID], data[genoID]
        if mygeno == '' : continue  # Skip blanks b/c were not found
        if myname in key and key[myname] != mygeno:
            print("\tWARNING! Line",myname,"given genotype",mygeno,"but already loaded with",key[myname])
        key[myname]=mygeno
    IN.close()

    print('\tLoaded genotype names for',len(key),"unique lines")
    return key

def update_blups(infile, key, outfile):
    print("Updating genotype names in", infile)
    IN = open(infile)
    OUT = open(outfile, "w")

    # Go through and change each name
    n, updated=0,0
    for line in IN:
        data=line.strip().split('\t')
        name=data[0]

        # Header lines; just write out unchanged
        if name in ['<Phenotype>', 'taxa', 'Taxon']:
            OUT.write(line)
            continue

        # Update name based on genotype key
        n+=1
        if name in key:
            data[0] = key[name]
            updated+=1
        else:   # If not in key, make sure there aren't any spaces or such that would cause problems
            print("\tNote: no genotype name for",name)
            data[0] = fix_bad_characters(name)
        OUT.write("\t".join(data)+'\n')

    print("\tUpdated",updated,"out of",n,"names")
    IN.close()
    OUT.close()

def fix_bad_characters(s):
    return s.replace(' ', '_')

if __name__ == '__main__': main()