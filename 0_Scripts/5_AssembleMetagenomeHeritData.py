__author__ = 'jgw87'

import argparse
import gzip
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import re

debug = False

def main():
    args = parse_args()
    if debug: args.infiles = args.infiles[:5]
    heritability, perms = calc_heritability(args.infiles)

    # Write out
    heritability.to_csv(args.outprefix + ".herits.txt", sep='\t')
    perms = pd.DataFrame(perms)
    perms.to_csv(args.outprefix + ".perms.txt", sep='\t')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", help="TASSEL MLM output files", nargs="*")
    parser.add_argument("-o", "--outprefix", help="Output text file prefix")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def calc_heritability(infiles):
    print("Calculating heritabilities in",len(infiles),"files")

    # First make a master dict() of all heritabilities; will split into perms and originals later
    heritability = dict()
    for infile in infiles:
        #print(infile)
        if len(open(infile).readline()) ==0 : continue #Skip empty files
        stats=pd.read_csv(infile, sep='\t')
        for trait, gen_var, resid_var in zip(stats['Trait'], stats['Genetic Var'], stats['Residual Var']):
            #if 'med1000_otu10_Rarefy10000_unifrac_bdiv_even10000_unweighted_unifrac_PC1' in trait:
                #print(trait,"has Vg", gen_var,"and Vr",resid_var)
            #print(trait)
            if trait in heritability: continue  # Already loaded, so skip
            heritability[trait] = gen_var / (gen_var + resid_var)

    # Now split into original traits and permutations
    core, perms = dict(), dict()
    for trait in sorted(heritability.keys()):
        if "_perm" in trait:
            basename = re.sub(string=trait, pattern="_perm.+", repl="")
            if basename not in perms:
                perms[basename] = list()
            perms[basename].append(heritability[trait])
        else:
            core[trait] = heritability[trait]


    traits = sorted(core.keys())
    for t in traits:
        if t not in perms: perms[t] = [np.nan, np.nan]  # Hack to get around traits with no permutation data
    h2 = [core[t] for t in traits]
    means = [np.nanmean(perms[t]) for t in traits]
    medians = [np.nanmedian(perms[t]) for t in traits]
    maxes = [np.nanmax(perms[t]) for t in traits]
    pvals = [np.sum(perms[t] >= core[t])/len(perms[t]) for t in traits] # Empirical p-value is how often got a heritability that high or more
    heritability = pd.DataFrame({"trait":traits, "h2":h2, "perm_mean":means, "perm_median":medians, "perm_max":maxes, "empirical_pval":pvals})
    heritability = heritability[["trait", "h2", "perm_mean", "perm_median","perm_max","empirical_pval"]]
    return heritability.sort(columns=["empirical_pval", "h2"], ascending=[True, False]), perms


if __name__ == "__main__":
    main()