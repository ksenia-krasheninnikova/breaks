#!/usr/bin/env python

import sys, os
from argparse import ArgumentParser
from breakpoint_graph import Anchor
from sumup_breaks import Breakpoint
import matplotlib
matplotlib.use('Agg')
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

DIVISOR = 250000000.

#colorkey from etc/colors.ucsc.conf
colors = {'scaffold0': '#483d8b',
'scaffold1'  : '#87ceeb',
'scaffold2'  : '#3cb371',
'scaffold3'  : '#7fff00',
#'scaffold4'  : '#ffff00',
'scaffold4'  : 'gold',
'scaffold5'  : '#daa520',
'scaffold6'  : '#cd5c5c',
'scaffold7' : '#f4a460',
'scaffold8'  : '#b22222',
'scaffold9' : '#ffa07a',
'scaffold10' : '#ff8c00',
'scaffold11' : '#ff4500',
'scaffold12' : '#ff1493',
'scaffold13' : '#b03060',
'scaffold14' : '#9370db',
'scaffold15' : '#20b2aa',
'scaffold16' : '#000080',
'scaffold17' : '#778899',
#'scaffold18' : '#ffe4c4',
'scaffold19' : '#556b2f'
}

canonical_names = {
'scaffold0' : 'A1', 
'scaffold1' : 'A2',
'scaffold2' : 'A3',
'scaffold3' : 'B1',
'scaffold4' : 'B2',
'scaffold5' : 'B3',
'scaffold6' : 'B4',
'scaffold7' : 'C1',
'scaffold8' : 'C2',
'scaffold9' : 'D1',
'scaffold10' : 'D2',
'scaffold11' : 'D3',
'scaffold12' : 'D4',
'scaffold13' : 'E1',
'scaffold14' : 'E2',
'scaffold15' : 'E3',
'scaffold16' : 'F1',
'scaffold17' : 'F2',
'scaffold19' : 'X'
}

def parse_sizes(path):
    sizes = {}
    with open(path) as f:
        for line in f:
            line = line.strip().split()
            if not line:
                return sizes
            sizes[line[0]] = int(line[-1])
    return sizes

def draw_chrom(ax, chrom_name, can_name, coords, x_start, width = 0.03):
    y_previous = 0
    rectangles = []
    prev_pos = 0
    for c in coords: 
        #print >> sys.stderr, chrom_name, c[0].pos, c[1].pos
        height = float(c[1].pos - c[0].pos) / DIVISOR
        #print >> sys.stderr,  colors[chrom_name]
        #if chrom_name == 'scaffold0' or chrom_name == 'scaffold16':
        #    r = Rectangle((x_start+0.03, y_previous + 0.04), width, height, facecolor = colors[chrom_name],label=can_name)
        #else:
        r = Rectangle((x_start+0.03, y_previous + 0.04), width, height, facecolor=colors[chrom_name], edgecolor='white', label='$\it{'+can_name+'}$')
        #print >> sys.stderr, 'height', c[1].pos - c[0].pos, 'diff', (c[0].pos-prev_pos)
        if not rectangles:
            rectangles.append(r)
        ax.add_patch(r)
        #y_previous += height + 0.03???
        y_previous += height
        prev_pos = c[1].pos
    return rectangles

def parse_breakpoints(path) :
    breakpoints = []
    with open(path) as f:
        for line in f:
            line = line.strip().split()
            if line[0] == 'None' or line[1] == 'None':
                continue
            breakpoints.append(Breakpoint(line[0],line[1],line[4]))
    return sorted(breakpoints, key=lambda x: x.left)

#make gap between breakpoints
def prepare_chrom(chrom, chrom_size, breaks, ref_hits):
    fbreaks = filter(lambda x: ref_hits[x.left].scaffold==chrom, breaks)
    left_breaks = map(lambda x: ref_hits[x.left], fbreaks)
    left_breaks = left_breaks + [Anchor(chrom, chrom_size, '')]
    right_breaks = map(lambda x: ref_hits[x.right], fbreaks)
    right_breaks = [Anchor(chrom, 1, '')] + right_breaks 
    return zip(right_breaks,left_breaks)
    

def run_drawings(breaks, sizes, ref_hits, suffix, name, verbose, threshold=0):
    pp=PdfPages('figure_chrom_'+suffix+'.pdf')
    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    width = 0.03
    x_start = 0.02
    rectangles = []
    names = []
    sorted_names = range(0,18) + [19] 
    sorted_names = map(lambda x: 'scaffold'+str(x), sorted_names)
    labels = []
    xticks = []
    for scaffold in sorted_names:
        xticks.append(x_start+0.03+float(width)/2)
        coords = prepare_chrom(scaffold, sizes[scaffold], breaks, ref_hits)
        if verbose:
            for c in coords:
                print c[0].to_string(), '-', c[1].to_string()
        #exit()
        can_name = canonical_names[scaffold]
        r = draw_chrom(ax, scaffold, can_name, coords, x_start)
        labels.append(can_name)
        rectangles += r
        #x_start += 1.5*width 
        x_start += 1.7*width 
    #leg=plt.legend(rectangles,labels,ncol=3, fontsize='xx-small', mode="expand", borderaxespad=0.)
    fig.patch.set_visible(False)
    #ax.axis('off')
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels)
    plt.tick_params(axis='x', labelsize=15, top='off', bottom='off')
    ax.axes.get_yaxis().set_visible(False)
    #fig.patch.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.title(' '.join(name.split('_')), fontsize=20)
    plt.tight_layout()
    plt.savefig(pp, format='pdf')    
    pp.close()
        
def load_hits(path):
    markers = {}
    with open(path) as f:
        for line in f:
            line = line.strip().split()
            markers[line[3][3:]] = Anchor(line[0], line[1], line[3][3:])
    return markers

#here we won't use filtered breakpoints without None on one of the side
#show breaks that should be done in reference in order to convert reference to the ${name} 
if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument('hits',nargs='?',default='./hits',help='directory with counted hits *counted.lst')
    parser.add_argument('breaks',help='breakpoints file')
    parser.add_argument('sizes',help='sizes of reference chromosomes')
    parser.add_argument('name', help='specie name to place in title')
    parser.add_argument('ref',nargs='?',default='FelCat')
    parser.add_argument('--threshold',nargs='?',default=0,type=int)
    parser.add_argument('-v',action='store_true')
    args = parser.parse_args()
    #print >> sys.stderr, 'loading markers...'
    #all_markers = load_markers(args.hits)
    print >> sys.stderr, 'loading reference hits...'
    ref_hits = load_hits(os.path.join(args.hits,'hits_'+args.ref+'.bed'))
    print >> sys.stderr, 'loading breakpoints...'
    breaks = parse_breakpoints(args.breaks)
    print >> sys.stderr, 'loading sizes...'
    sizes = parse_sizes(args.sizes)
    print >> sys.stderr, 'drawings...'
    suffix = '' if not args.threshold else '_'+str(args.threshold)
    suffix = args.name + suffix
    run_drawings(breaks, sizes, ref_hits, suffix, args.name, args.v, args.threshold)
    
