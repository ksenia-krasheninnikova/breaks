#!/usr/bin/env python

from argparse import ArgumentParser

from pycbio.hgdata.psl import PslReader


def convert_blocks(psl_file):
    queries = []
    targets = []
    for psl in PslReader(psl_file):
        query = '\t'.join(map(str, [psl.qName, psl.qStart, psl.qEnd]))
        queries.append(query)
        target = '\t'.join(map(str, [psl.tName, psl.tStart, psl.tEnd]))
        targets.append(target)
    return queries, targets

def print_out(blocks, filename):
    with open(filename, 'w') as f:
        for b in blocks:
            f.write(b+'\n')

def print_out_together(psl_file):
    for psl in PslReader(psl_file):
        print '\t'.join(map(str, [psl.qName, psl.qStart, psl.qEnd, \
            psl.tName, psl.tStart, psl.tEnd]))

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('psl')
    parser.add_argument('queries',nargs='?', help='bed file for queries')
    parser.add_argument('targets',nargs='?', help='bed file for targets')
    parser.add_argument('--one_file', action='store_true', help='put query and target in one file')
    args = parser.parse_args()
    if args.one_file:
        print_out_together(args.psl)
    else:
        queries,targets = convert_blocks(args.psl)
        print_out(queries, args.queries)
        print_out(targets, args.targets)
