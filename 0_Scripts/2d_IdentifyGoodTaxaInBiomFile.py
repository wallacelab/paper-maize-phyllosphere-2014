__author__ = 'jgwall'

import argparse

debug = False


def main():
    args = parse_args()
    print("Finding good taxa from biom summary file",args.infile)

    samples, bad_taxa = dict(), set(args.remove_samples)

    # Load in all the samples and their depth
    IN = open(args.infile, "r")
    skipped=0
    while(IN.readline().strip() != "Counts/sample detail:"):    # Proceed through until hit relevant part of file
        pass
    for line in IN:
        sample, depth = line.strip().split(" ")
        sample = sample.rstrip(":") # Remove colon
        if sample in bad_taxa:
            skipped+=1
            continue
        samples[sample] = float(depth)
    print("\tLoaded",len(samples),"samples. (Skipped",skipped," because flagged as bad from the command line)")
    if skipped != len(bad_taxa):
        print("\t\tWARNING! Not all taxa flagged as bad could be found! Expected",len(bad_taxa),"but found",skipped);
    print("\t\tExpected bad samples:",args.remove_samples);

    # Find the highest read depth among blank samples
    blank_depths = list()
    for s in samples:
        for pattern in args.blank_patterns:
            if pattern in s:
                blank_depths.append(samples[s])
                break
    max_blank = max(blank_depths)
    print("\tFound",len(blank_depths),"blanks with a maxmimum depth of",max_blank)

    # Determine bad samples based off read depth
    n=0
    for s in samples:
        if samples[s] < max_blank * args.sample_to_blank_threshhold:
            n+=1
            bad_taxa.add(s)
    print("\tFound",n,"taxa with too little depth to pass. (This includes the blanks themselves)")


    if(args.outfile):
        OUT = open(args.outfile, "w")
        n=0
        for s in sorted(samples.keys()):
            if s not in bad_taxa:
                n+=1
                OUT.write(s + "\n")
        OUT.close()
        print("Total of",n,"good taxa identified. Output to",args.outfile)


    if(args.failure_file):
        print("Total of",len(bad_taxa),"bad taxa identified. Outputting to",args.failure_file)
        OUT = open(args.failure_file, "w")
        for s in sorted(bad_taxa):
            OUT.write(s + "\n")
        OUT.close()



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="biom summary file (across samples, not observations)")
    parser.add_argument("-o", "--outfile", help="Output file to list the names of the good taxa")
    parser.add_argument("-f", "--failure-file", help="Output file to list the names of the bad taxa")
    parser.add_argument("-b", "--blank-patterns", nargs="*", help="Key words to recognize blanks in the sample names")
    parser.add_argument("-t", "--sample-to-blank-threshhold", type=float, default=0.0,
                        help="How many more reads a sample has to have over the most abundant blank to pass filtering")
    parser.add_argument("-r", "--remove-samples", nargs="*", help="List of samples that should be removed _before_ looking for blanks, etc.")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


if __name__ == '__main__': main()