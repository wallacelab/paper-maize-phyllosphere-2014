__author__ = 'jgwall'

import argparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import re

debug = False
matplotlib.rcParams.update({'font.size': 10})

def main():
    args = parse_args()
    print("Plotting alpha diversity")
    data = [load_averages(i) for i in args.infiles]

    # Collate by metric and category
    metrics = sorted({d['metric'] for d in data})
    categories = sorted({d['category'] for d in data})

    # Make plots
    nrow = len(metrics)
    ncol=len(categories)
    fig = plt.figure(figsize=(4*ncol, 3*nrow))
    grid =gridspec.GridSpec(nrows=nrow, ncols=ncol, hspace=0.25, wspace=0.25)

    # Loop over datasets and plot out
    plotdata={'metric', 'category', 'xticks', 'xmax'}   # Which keys hold axis data
    for d in data:
        myrow, mycol = metrics.index(d['metric']), categories.index(d['category'])
        ax = fig.add_subplot(grid[myrow, mycol])

        # Add individual datasets
        for group in sorted(d.keys()):
            if group in plotdata: continue  # Skip plot metadata
            mydata=d[group]
            ax.errorbar(x=d['xticks'], y=mydata['points'], yerr=mydata['error'], label=group, color=mydata['color'])

        # Integrate metadata
        ax.set_xlim([0, d['xmax']])
        ax.set_title(d['category'], weight='bold')
        if(mycol==0): ax.set_ylabel(d['metric'], weight='bold') # Only set y label on leftmost plot

        # Make pretty legend
        handles, labels = ax.get_legend_handles_labels()
        handles = [h[0] for h in handles]   # Remove errorbars
        ax.legend(handles, labels, fontsize='xx-small', loc='upper left', frameon=False, ncol=2 if len(d)>10 else 1)

    fig.savefig(args.outfile, dpi=150)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", nargs="*")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def load_averages(infile):
    print("\tLoading average alpha diversity data from",infile)
    IN = open(infile)

    # Get metric and category
    metric=IN.readline().lstrip("#").strip()
    metric = re.sub(string=metric, pattern="\.txt$", repl="")
    category = IN.readline().lstrip("#").strip()

    # Get axis values
    xaxis = IN.readline().strip().split()[1:]
    xaxis = [float(x) for x in xaxis]

    # Get x max
    xmax =  IN.readline().strip().split()
    xmax = float(xmax[1])

    # Store information for later
    data = dict()
    data['metric']=metric
    data['category']=category
    data['xticks'] = xaxis
    data['xmax'] = xmax

    # Parse individual datasets
    for name in IN:
        color = IN.readline().strip().split()[1]
        series = IN.readline().strip().split()[1:]
        error = IN.readline().strip().split()[1:]

        # Clean name
        name = name.lstrip(">").strip()
        if name == '8082014': name='8 Aug 2014'
        if name == '8262014': name = '26 Aug 2014'
        if name == 'nss': name='Non-stiff-stalk'
        if name == 'ss': name = 'Stiff-stalk'
        if name == 'ts': name = 'Tropical'
        name=name.title()

        # Convert plot values to floats
        series = [float(s) for s in series]
        error = [float(e) for e in error]

        data[name]=dict()
        data[name]['color']=color
        data[name]['points'] = series
        data[name]['error'] = error
    print('\t\tLoaded plotting data for',len(data),"items")
    return(data)


if __name__ == '__main__': main()