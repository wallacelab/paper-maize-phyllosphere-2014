__author__ = 'jgwall'

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re


debug = False
taxa_delim="|"    # character to delimit levels of taxonomy

def main():
    args = parse_args()
    print("Plotting OTU distributions within and across samples")

    # Load data
    biom=pd.read_table(args.infile, index_col=0)
    biom.index = [str(i) for i in biom.index]   # Convert it all to strings to make downstream stuff easier
    taxonomy = load_taxonomy(args.taxonomyfile)
    colorkey=pd.read_table(args.colorkey)
    core_otus=[str(o) for o in pd.read_table(args.core_otus)['otu']]
    print(core_otus)

    # Get per-sample sums based on the requested taxa divisions
    divisions = divide_otus_by_taxon(colorkey['taxon'], biom, taxonomy)
    core = np.sum(biom.loc[core_otus,:], axis=0)

    # Convert into fractions of total
    totals = np.sum(biom, axis=0)
    core /= totals
    for taxon in divisions: divisions[taxon] /= totals

    # Write out core fraction to text
    outtext=args.outprefix + ".core_fractions.txt"
    print("Writing out fraction core OTUs to", outtext)
    core.to_csv(outtext, sep='\t')
    print("\tCore OTU stats:")
    print("\t\tmean:",np.mean(core))
    print("\t\tmedian:", np.median(core))
    print("\t\tmin:", np.min(core))
    print("\t\tmax:", np.max(core))


    # Plot graphic
    fig = plt.figure(figsize=(12,6))

    ax_across = fig.add_axes((0.2, 0.2, 0.3, 0.75))
    plot_dist_across_samples(ax_across, core, divisions, colorkey)

    ax_within = fig.add_axes((0.53, 0.2, 0.4, 0.75))
    plot_bars_within_samples(ax_within, divisions, colorkey, args.num_samples, args.sample_filter, args.seed)

    fig.savefig(args.outprefix + ".png", dpi=100)
    fig.savefig(args.outprefix + ".svg", dpi=100)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Text OTU table of abundances")
    parser.add_argument("-t", "--taxonomyfile", help="Uclust-assigned taxonomy")
    parser.add_argument("-c", "--colorkey", help="Color key for various taxa. NOTE: Order is important, as higher-level taxa will exlclude those included in prior lower-level taxa")
    parser.add_argument("--core-otus", help="File listing core OTUs")
    parser.add_argument("-o", "--outprefix")
    parser.add_argument("-s", "--seed", type=int, default=1, help="Random seed for choosing samples for barplot")
    parser.add_argument("-n", "--num-samples", type=int, default=10, help="Number of random samples to plot for barplot")
    parser.add_argument("--sample-filter", help="Pattern to filter samples by before subsetting. (Only ones with this pattern are kept.)")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def load_taxonomy(infile):
    print("\tLoading taxonomy")
    data = pd.read_csv(infile, sep='\t', index_col=0, header=None, dtype=object)
    # print(data.head())
    data.columns = ["taxonomy","unknown1","unknown2"]
    taxonomy=dict()
    for i in range(len(data)):
        otu = str(data.index[i])
        taxonomy[otu] = taxonomy_from_string(data["taxonomy"].iloc[i], otu)
        taxonomy[otu] = taxa_delim.join(taxonomy[otu])    # turn into a single string
        # print(otu,"has taxonomy",taxonomy[otu])
    return taxonomy

def taxonomy_from_string(string, otu):
    notfound=""
    kingdom = re.search(string=string, pattern="k__([^;]*)").group(1) if "k__" in string else notfound
    phylum = re.search(string=string, pattern="p__([^;]*)").group(1) if "p__" in string else notfound
    class_ = re.search(string=string, pattern="c__([^;]*)").group(1) if "c__" in string else notfound
    order = re.search(string=string, pattern="o__([^;]*)").group(1) if "o__" in string else notfound
    family = re.search(string=string, pattern="f__([^;]*)").group(1) if "f__" in string else notfound
    genus = re.search(string=string, pattern="g__([^;]*)").group(1) if "g__" in string else notfound
    species = re.search(string=string, pattern="s__([^;]*)").group(1) if "s__" in string else notfound

    formatted=[kingdom]
    if phylum != notfound: formatted.append(phylum)
    if class_ != notfound: formatted.append(class_)
    if order != notfound:  formatted.append(order)
    if family != notfound: formatted.append(family)
    if genus != notfound:  formatted.append(genus)
    if species != notfound: formatted.append(species)
    formatted.append(otu)

    # print(string)
    # print(formatted,"\n")

    return formatted

def divide_otus_by_taxon(taxa, biom, taxonomy):
    print("\tCalculating # of reads for taxonomic divisions")
    divisions=dict()
    loaded=set()    # Keep track of which OTUs are already loaded
    for taxon in taxa:
        all_targets = [taxon in taxonomy[str(otu)] for otu in biom.index]
        new_targets = np.array(all_targets) & np.array([otu not in loaded for otu in biom.index])
        for o in biom.index[new_targets]:
            loaded.add(o)

        print("\t\t",taxon,"found for",np.sum(all_targets),"OTUs, of which",np.sum(new_targets),"are not already assigned to a group")

        # Subset out and take the per-sample sum of OTUs
        subbiom = biom.loc[new_targets, :]
        divisions[taxon] = np.sum(subbiom, axis=0)

    # Sanity check that numbers all add up
    orig_totals = np.sum(biom, axis=0)
    reconstructed=pd.DataFrame(divisions)
    new_totals = np.sum(reconstructed, axis=1)
    mismatches=0
    for old, new in zip(orig_totals, new_totals):
        if old!=new: mismatches+=1
        # print(old, new)
    if mismatches>0:
        print("\t\t### WARNING! ### Found",mismatches,"instances of new totals not matching the ones from the original biom file")
        print("\t\t\tThis may be because the target taxa are not exhaustive; best to check your file to be sure")

    return divisions

# barplots of distributions within samples
def plot_bars_within_samples(ax, divisions, colorkey, num_samples, sample_filter, seed):
    print("Plotting barplots of samples")

    # Make list of all samples
    samples=set()
    for d in divisions:
        for s in divisions[d].index:
            samples.add(s)
    print('\t',len(samples),"samples found to select for printing")

    # Filter out if required
    if sample_filter:
        toremove=set()
        for s in samples:
            if sample_filter not in s:
                toremove.add(s)
        samples-=toremove
        print("\tAfter filtering for pattern '"+sample_filter+"', have",len(samples),"samples left to choose from")

    # Select random subset and filter
    np.random.seed(seed)
    samples=sorted(samples)
    target_samples = sorted(np.random.choice(samples, size=num_samples, replace=False))
    bardata = dict()
    for group in divisions:
        bardata[group] = divisions[group][target_samples]

    # Make graphics
    bottom = np.zeros(num_samples)
    xvals = range(num_samples)
    bar_width=0.7
    for taxon in colorkey['taxon']:
        myvals = bardata[taxon]
        mycolor = colorkey['color'][colorkey['taxon'] == taxon]
        mylabel = colorkey['display_name'][colorkey['taxon'] == taxon].iloc[0]
        ax.bar(left=xvals, height=myvals, width=bar_width, bottom=bottom, color=mycolor, align='center', label=mylabel, linewidth=0.5)
        bottom += myvals

    # # Add legend
    legend_properties = {'weight': 'bold', 'size': 'xx-small'}
    ax.legend(prop=legend_properties, ncol=2, frameon=False, loc='upper center')

    xticklabels = target_samples
    if sample_filter:
        xticklabels = [re.sub(pattern=sample_filter, string=x, repl="") for x in xticklabels]

    # Prettify x axis
    xvals = list(xvals)
    xpad=1
    ax.set_xlim([min(xvals)-xpad, max(xvals)+xpad])
    ax.set_xticks(xvals)
    ax.set_xticklabels(xticklabels, weight='bold', size='small', rotation='vertical')
    # ax.set_xticklabels([" "] * len(samples))
    ax.xaxis.set_ticks_position('none')  # bottom, top, both, or none

    # Prettify y axis
    yticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    ax.set_ylim(-0.05, 1.25)
    ax.yaxis.set_ticks_position('right')		# left, right, both, or none
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(y) for y in yticks], weight='bold', size='small')
    ax.yaxis.set_label_position("right")
    ax.set_ylabel("Fraction of reads in each sample", weight='bold', labelpad=16, rotation=-90)
    ax.set_xlabel("Individual Samples", weight='bold')


# Violinplots of distribtuion of taxa
def plot_dist_across_samples(ax, core, divisions, colorkey):
    print("Plotting violint plots of distributions across samples")
    violins=dict()
    xvals = np.arange(len(colorkey)) #+ 1

    # Plot core
    violins['core'] = ax.violinplot(np.array(core), positions=[max(xvals) + 1], showextrema=False, vert=False)

    # Plot others
    for x, taxon in zip(reversed(xvals), reversed(colorkey['taxon'])):
        violins[taxon] = ax.violinplot(divisions[taxon], positions=[x], showextrema=False, vert=False)

    # Colorize violin plots
    for taxon in violins:
        bodies = violins[taxon]['bodies']

        if taxon == "core":
            mycolor='black'
            mylabel="Core OTUs"
        else:
            mycolor = colorkey['color'][colorkey['taxon']==taxon]
            mylabel = colorkey['display_name'][colorkey['taxon'] == taxon].iloc[0]
        for b in bodies:
            b.set_color(mycolor)
            # b.set_color('black')
            b.set_alpha(1)
            # b.set_linewidth(0.5)
            b.set_edgecolor('black')
            b.set_linewidth(0.1)
            b.set_label(mylabel)

    # # Add legend
    # legend_properties = {'weight': 'bold', 'size': 'xx-small'}
    # ax.legend(prop=legend_properties, ncol=2, frameon=False, mode='expand')

    # Prettify axes
    yvals = list(xvals) + [max(xvals)+1]
    xticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0]
    ax.set_yticks(yvals)
    ax.set_yticklabels(list(colorkey['display_name']) + ["Core OTUs"], weight='bold', size='small')
    # ax.set_yticklabels([""] * len(yvals))
    ax.set_xlim(-0.05, 1)
    # ax.set_ylim([min(xvals)-1, max(xvals)+4])
    ax.yaxis.set_ticks_position('none')	# bottom, top, both, or none
    ax.xaxis.set_ticks_position('bottom')		# left, right, both, or none
    ax.set_xlabel("Fraction of reads in each sample", weight='bold')
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(y) for y in xticks], weight='bold', size='small')

# Older, vertical version
# def plot_dist_across_samples(ax, core, divisions, colorkey):
#     print("Plotting violint plots of distributions across samples")
#     violins=dict()
#
#     # Plot core
#     violins['core'] = ax.violinplot(np.array(core), positions=[0], showextrema=False)
#
#     # Plot others
#     xvals = np.arange(len(colorkey)) + 1
#     for x, taxon in zip(xvals, colorkey['taxon']):
#         violins[taxon] = ax.violinplot(divisions[taxon], positions=[x], showextrema=False)
#
#     # Colorize violin plots
#     for taxon in violins:
#         bodies = violins[taxon]['bodies']
#         if taxon == "core":
#             mycolor='black'
#         else:
#             mycolor = colorkey['color'][colorkey['taxon']==taxon]
#         for b in bodies:
#             b.set_color(mycolor)
#             # b.set_color('black')
#             b.set_alpha(1)
#             b.set_linewidth(0)
#
#     # Prettify axes
#     xvals = [0] + list(xvals)
#     yticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0]
#     ax.set_xticks(xvals)
#     ax.set_xticklabels(["Core OTUs"] + list(colorkey['display_name']), rotation='vertical', weight='bold', size='small')
#     ax.set_ylim(-0.05, 1)
#     ax.xaxis.set_ticks_position('bottom')	# bottom, top, both, or none
#     ax.yaxis.set_ticks_position('left')		# left, right, both, or none
#     ax.set_ylabel("Fraction of reads in each sample", weight='bold')
#     ax.set_yticks(yticks)
#     ax.set_yticklabels([str(y) for y in yticks], weight='bold', size='small')




if __name__ == '__main__': main()