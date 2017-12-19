#! /usr/bin/python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--source")
parser.add_argument("-d", "--destination")
parser.add_argument("-o", "--outfile")
args=parser.parse_args()

def file_to_dict(infile):
  print("Loading data from",infile)
  data = open(infile, 'r').readlines()
  data = [line.strip().split('\t') for line in data]
  key = { d[0]: d[1] for d in data}
  print("\tResulting dictionary of MD5sums has",len(key),"items")
  if(len(key) != len(data)):
    print("Warning! Length of input was",len(data),"; some items may have been lost due to duplicated names!")
  return key

source = file_to_dict(args.source)
dest = file_to_dict(args.destination)

# Check for mismatches
mismatches=set()
for file in source:
  if file not in dest:	# Skip ones that don't match; will handle these later
    print(file)
    continue
  if source[file] != dest[file]:
    mismatches.add(file)

# Check for missing files
source_set = set(source.keys())
dest_set = set(dest.keys())
dest_lacking = source_set.difference(dest_set)
source_lacking = dest_set.difference(source_set)

# Print report
print("Report:")
print("\t",len(dest_lacking),"files are in the source location but not the destination:",sorted(dest_lacking))
print("\t",len(source_lacking),"files are in the destination location but not the source:",sorted(source_lacking))
print("\t",len(mismatches),"files have different MD5sums in the two locations:",sorted(mismatches))

# Output mistmatches
OUT = open(args.outfile, "w")
OUT.write("file\tsource\tdest\n")
for name in sorted(mismatches):
  OUT.write("\t".join([name, source[name], dest[name]]) + "\n")
OUT.close()