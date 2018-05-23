#!/usr/bin/env python

from argparse import ArgumentParser
from itertools import groupby
import time

from pycbio.hgdata.psl import PslReader
from pycbio.hgdata.psl import Psl

def get_blocks_set(psl_file):
    blocks = []
    for psl in PslReader(psl_file):
        blocks += psl.blocks
    return set(blocks)

#assume a.start < b.start
def are_syntenic(a, b):
    return a.qEnd <= b.qStart and\
            a.tEnd <= b.tStart and\
             a.psl.tName == b.psl.tName and\
                a.psl.strand == b.psl.strand


def is_not_overlapping_ordered_pair(a, b, threshold=5000):
    return are_syntenic(a, b) and\
             0 <= b.qStart - a.qEnd < threshold and\
              0 <=  b.tStart - a.tEnd < threshold


def get_next(pos, query_group, max_anchor_distance=5000):
    f = []
    for i in range(pos+1, len(query_group)):
        if is_not_overlapping_ordered_pair(query_group[pos], query_group[i], max_anchor_distance):
            if not f:
                f.append(i)
            elif f and is_not_overlapping_ordered_pair(query_group[f[0]], query_group[i], max_anchor_distance):
                return f
            else:
                f.append(i)
    return f


'''
dag is a dict that for a given vertex stores all its possible next verties
hidden_vertices is a set of vertices that are already in paths
'''
def weigh_dag(group, dag, hidden_vertices, max_anchor_distance):
    #weight of the edge equals length 
    #of the next block
    #weight of a vertice equals 
    #estimated weight: w_j < w_i + w_e(ij) =>
    #updated w_j
    #also remember how we came here: 
    #(prev_vertex, weight)
    weighted_dag = {}
    for i in range(len(group)):
        if i in hidden_vertices:
            continue
        #print i, group[i], len(group)
        if not i in dag:
            nexts = get_next(i, group, max_anchor_distance)
            dag[i] = nexts
        else:
            nexts = dag[i]
        #if never visited this vertex then 
        #its weight equals to its size
        #because otherwise we will never count its size
        if not i in weighted_dag:
            weighted_dag[i] = (-1, group[i].size)
        for j in nexts:
            if j in hidden_vertices:
                continue
            alternative_weight = weighted_dag[i][1] + group[j].size
            if not j in weighted_dag or weighted_dag[j][1] < alternative_weight: 
                #w_i + weight of the next edge 
                weighted_dag[j] = (i, alternative_weight)
    return weighted_dag

def traceback(weighted_dag, hidden, group):
    #get the heaviest path weight
    start_vertex = max(weighted_dag.items(), key=lambda x:x[1][1])[0]
    path = [start_vertex]
    prev_vertex = weighted_dag[start_vertex][0]
    while prev_vertex != -1:
        path.append(prev_vertex)
        prev_vertex = weighted_dag[prev_vertex][0]
    hidden.update(set(path))
    return map(lambda x: group[x], path)[::-1]

def dag_merge(psl, min_block_breath, max_anchor_distance):
    blocks = get_blocks_set(psl)
    blocks = sorted(blocks, key=lambda x: x.psl.qName)
    blocks_grouped_by_query = map(lambda x:list(x[1]), groupby(blocks, key=lambda x:x.psl.qName))
    paths = []
    for group in blocks_grouped_by_query:
        dag = {}
        hidden_vertices = set()
        group = sorted(group, key=lambda x: x.qStart)
        while len(group) != len(hidden_vertices):
            ts = time.time()
            weighted_dag = weigh_dag(group, dag, hidden_vertices, max_anchor_distance)
            path = traceback(weighted_dag, hidden_vertices, group)
            if not path:
                break
            qLen = path[-1].qEnd - path[0].qStart
            tLen = path[-1].tEnd - path[0].tStart
            if qLen >= min_block_breath and tLen >= min_block_breath:
                paths.append(path)
            #think if we should keep path that is not complies with lengths bounds
            #we could traverse these vertices later in other paths
            #for e in path:
            #    group.remove(e)
    return paths
            
def construct_psl(blocks):
    psl = Psl()
    psl.match = sum(map(lambda x: x.qEnd-x.qStart, blocks))
    psl.misMatch = 0
    psl.repMatch = 0
    psl.nCount = 0
    psl.qNumInsert = len(filter(lambda x: x>0, map(lambda x: x[1].qStart-x[0].qEnd, zip(blocks,blocks[1:]))))
    psl.qBaseInsert = sum(map(lambda x: x[1].qStart-x[0].qEnd, zip(blocks,blocks[1:])))
    psl.tNumInsert = len(filter(lambda x: x>0, map(lambda x: x[1].tStart-x[0].tEnd, zip(blocks,blocks[1:]))))
    psl.tBaseInsert = sum(map(lambda x: x[1].tStart-x[0].tEnd, zip(blocks,blocks[1:])))
    psl.qName = blocks[0].psl.qName
    psl.qSize = blocks[0].psl.qSize
    psl.qStart = blocks[0].qStart
    psl.qEnd = blocks[-1].qEnd
    psl.tName = blocks[0].psl.tName
    psl.tSize = blocks[0].psl.tSize
    psl.strand = blocks[0].psl.strand
    if psl.strand == '++':
        psl.tStart = blocks[0].tStart
        psl.tEnd = blocks[-1].tEnd
    elif psl.strand == '+-':
        psl.tEnd = psl.tSize - blocks[0].tStart
        psl.tStart = psl.tSize - blocks[-1].tEnd
    psl.blockCount = len(blocks)
    for b in blocks:
        b.psl = psl
    psl.blocks = blocks
    return psl
    
if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('min_block_size', nargs='?', default=5000, type=int)
    parser.add_argument('max_anchor_distance', nargs='?', default=5000, type=int)
    parser.add_argument('psl')
    parser.add_argument('out')
    args = parser.parse_args()
    print 'dag merge...'
    print 'using min_block_size = ', args.min_block_size, \
        'max_anchor_distance =', args.max_anchor_distance
    merged = dag_merge(args.psl, args.min_block_size, args.max_anchor_distance)
    print 'storing output...'
    with open(args.out, 'w') as f:
        for blocks in merged:
           psl = construct_psl(blocks)
           f.write('\t'.join(psl.toRow())+'\n')
