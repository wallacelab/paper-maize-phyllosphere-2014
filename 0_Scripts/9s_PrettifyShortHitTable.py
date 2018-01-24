__author__ = 'jgwall'

import argparse
import pandas as pd
import re

debug = False


def main():
    args = parse_args()
    print("Prettifying short hit table in",args.infile)

    # Load metagenome key
    raw_key = pd.read_table(args.metagenome_key)
    key=dict()
    for id, desc in zip(raw_key['id'], raw_key['description']):
        key[id]=desc

    # Parse table header
    IN = open(args.infile, "r")
    header=IN.readline()

    # Parse table contents
    n=0
    output=list()
    for line in IN:
        n+=1
        # Turn into a nested list
        mydata = line.strip().split("\t")
        mydata = [d.split(',') for d in mydata]
        # print(mydata)

        # Prettify and sort traits
        for i in range(len(mydata)):
            for j in range(len(mydata[i])):
                mydata[i][j] = prettify_name(mydata[i][j], key)
            mydata[i] = sorted(mydata[i])   # Sort

        # Now go through and print in order, substituting blanks if a sub-list is not long enough for the longest sub-list
        lengths = [len(d) for d in mydata]
        maxlength = max(lengths)
        for i in range(maxlength):
            myline = list()
            # Walk through the columns on this line
            for j in range(len(mydata)):
                if len(mydata[j]) <= i:
                    myline.append("")
                else:
                    myline.append(mydata[j][i])
            output.append("\t".join(myline) + "\n")
    print("\tConverted",n,"input lines into",len(output),"output lines")

    # write out
    output = [header] + output  # Add header to beginning
    OUT = open(args.outfile, "w")
    OUT.writelines(output)
    OUT.close()

# Helper function to clean up names for table; includes general rules and some specific cases
def prettify_name(x, key):
    # First, check for pattern that matches marker names and skip altering them any.
    if re.match(string=x, pattern="^S[0-9]+_[0-9]+$"): return x

    x = re.sub(string=x, pattern="_[0-9]+$", repl="")

    # Add metagenome info
    name = re.sub(string=x, pattern="^log_", repl="")
    if name in key:
        description = key[name]
        x = x + ": " + cap_first(description)

    # specific cases
    x = re.sub(string=x, pattern="_", repl=" ")
    x = re.sub(string=x, pattern="Meiosis   yeast", repl="Meiosis - yeast")
    x = re.sub(string=x, pattern="D Arginine and D ornithine", repl="D-Arginine and D-ornithine")
    x = re.sub(string=x, pattern="qiime.+Rhizobiales", repl="Rhizobiales (order)")
    x = re.sub(string=x, pattern="unifrac", repl="UniFrac")

    # Move log-transofrmed to end
    x = re.sub(string=x, pattern="^log (.+)", repl="\\1 (log-transformed)")

    # Capitalize the first letter
    x = cap_first(x)

    return(x)

# Helper function to caoptialize the first letter of something
def cap_first(x):
    return x[0].upper() +  x[1:]

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-m", "--metagenome-key")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


if __name__ == '__main__': main()