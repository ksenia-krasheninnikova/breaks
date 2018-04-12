#!/usr/bin/env python

from argparse import ArgumentParser
from collections import defaultdict


from model import Anchor, empty_anchor, parse_counted

def parse_seqs(path):
    seqs = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            seqs.add(line)
    return seqs


def parse_anchors(path, seqs=False, list_of_IDs=set()):
    if seqs:
        seqs = parse_seqs(seqs)
    genome = defaultdict(list)
    with open(path) as f:
        for line in f:
            line = line.strip().split()
            if seqs and line[0] not in seqs:
                continue
            ID = line[3][3:]
            if list_of_IDs and (not ID in list_of_IDs):
                continue
            genome[line[0]].append(Anchor(line[0], int(line[1]), ID))
    for chrom in genome:
        genome[chrom] = sorted(genome[chrom],key=lambda x: x.pos)
        genome[chrom].insert(0,empty_anchor)
        genome[chrom].append(empty_anchor)
    return genome


def report_breakpoint(common_anc, qanc, rancid):
    print(str(common_anc) + ' - ' + str(qanc) + ' | ' +\
        common_anc.ID + ' - ' + rancid)


def print_genome(g):
    for scaffold in g:
        print scaffold
        for anchor in g[scaffold]:
            print anchor.ID,
        print


def get_query_pos(ID, q):
    for l in q:
        f = filter(lambda x: x.ID == ID, l) 
        if f:
            ind = l.index(f[0])
            return l, ind
    raise Exception('No such query ID '+ID)


def make_index(r):
    ind = {}
    for scaffold in r:
        for anc in r[scaffold]:
            ind[anc.ID] = anc
    return ind


def print_r_as_q(r_ind, q):
    for scaffold in q:
        print scaffold
        for anc in q[scaffold]:
            if not anc.ID:
                continue
            print anc.ID, r_ind[anc.ID].scaffold, r_ind[anc.ID].pos
    

def print_q_as_r(r, q):
    for scaffold in r:
        for ref_anc in r[scaffold]:
            if not ref_anc.ID:
                continue
            l, ind = get_query_pos(ref_anc.ID, q)
            query_anc = l[ind]
            position = 'Middle'
            if not l[ind+1].ID:
                position = 'Last'
            if not l[ind-1].ID:
                position = 'First'
            if len(l) == 1:
                position = 'Only'
            print(str(ref_anc) + ' - ' + str(query_anc) + ' ' + position)


def run(r, q):
    for scaffold in r:
        for i in range(1,len(r[scaffold])-2):
            ref_triple = (r[scaffold][i-1],r[scaffold][i],r[scaffold][i+1])
            l, ind = get_query_pos(r[scaffold][i].ID, q)
            query_triple = (l[ind-1],l[ind],l[ind+1])
            ref_neighb = [ref_triple[0].ID, ref_triple[2].ID]
            left_break = right_break = False
            if query_triple[0].ID != None and not (query_triple[0].ID in ref_neighb):
                if query_triple[2].ID:
                    #collect the breakpoint in reference too
                    if query_triple[2].ID in ref_neighb:
                        ref_br = filter(lambda x: x!=query_triple[2].ID, ref_neighb)
                        ref_br = ref_br[0]
                    else:
                        ref_br = ref_neighb[0]
                #if one of neighbours in reference is not None
                elif bool(ref_neighb[0] == None) != bool(ref_neighb[1] == None) :
                    ref_br = filter(lambda x: x != None, ref_neighb)
                    ref_br = ref_br[0]
                else:
                    ref_br = None
                report_breakpoint(query_triple[1], query_triple[0], str(ref_br))
                left_break = True
            if query_triple[2].ID != None and not (query_triple[2].ID in ref_neighb):
                if query_triple[0].ID:
                    if query_triple[0].ID in ref_neighb:
                        ref_br = filter(lambda x: x!=query_triple[0].ID, ref_neighb)
                        ref_br = ref_br[0]
                    else:
                        ref_br = ref_neighb[1]
                elif bool(ref_neighb[0] == None) != bool(ref_neighb[1] == None) :
                    ref_br = filter(lambda x: x != None, ref_neighb)
                    ref_br = ref_br[0]
                else:
                    ref_br = None
                report_breakpoint(query_triple[1], query_triple[2], str(ref_br))
                right_break = True
            if left_break and right_break:
                new_l = [l[:ind]+[empty_anchor], [empty_anchor, l[ind], empty_anchor], [empty_anchor]+l[ind+1:]]
            elif not left_break and right_break:
                new_l = [l[:ind+1]+[empty_anchor], [empty_anchor]+l[ind+1:]]
            elif left_break and not right_break:
                new_l = [l[:ind]+[empty_anchor], [empty_anchor]+l[ind:]]
            else:
                new_l = [l]
            q.remove(l)
            q += new_l 


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--reference',help='file with reference anchors')
    parser.add_argument('--query',help='file with liftovered anchors')
    parser.add_argument('--print_ref',action='store_true')
    parser.add_argument('--print_q',action='store_true')
    parser.add_argument('--print_q_as_r',action='store_true',help='print reference constructed from query')
    parser.add_argument('--print_r_as_q',action='store_true',help='print query constructed from reference')
    parser.add_argument('--seqs',help='pass sequences names on query to look for breakpoints')
    args = parser.parse_args()
    counted = parse_counted(args.query)
    r = parse_anchors(args.reference, seqs=args.seqs, list_of_IDs=counted)
    ind = make_index(r)
    q = parse_anchors(args.query, list_of_IDs=set(ind.keys()))
    if args.print_q:
        print_genome(q)
        exit()
    if args.print_r_as_q:
        print_r_as_q(ind, q)
        exit()
    q = q.values()
    if args.print_ref:
        print_genome(r)
        exit()
    if args.print_q_as_r:
        print_q_as_r(r, q)
        exit()
    run(r, q)

