#!/usr/bin/env python

import os
from argparse import ArgumentParser

from cluster_utils import run_joblist

prefix='/hive/groups/recon/projs/felidae_comp/synteny-play-ground'

def make_joblist(specie, common_name, max_jobs=4):
    chromosomes = range(0,20)
    chromosomes = map(str, chromosomes)
    #chromosomes = map(lambda x: 'scaffold'+x, chromosomes)
    joblist=os.path.join(prefix,'data/jobs.txt')
    with open(joblist, 'w') as f:
       for c in chromosomes:
            script = os.path.join(prefix,'bin/psl_merger_wrapper.sh')
            folder=os.path.join(prefix,'data/felidae/',common_name)
            if not os.path.isdir(folder):
                raise Exception('No input data folder', folder)
            data = os.path.join(folder,specie+'.FelisCatus.'+c+'.psl')
            folder=os.path.join(prefix,'data/human/',common_name)
            if not os.path.isdir(folder):
                os.makedirs(folder)
            output = os.path.join(folder,specie+'.FelisCatus.'+c+'.merged.psl')
            line = ' '.join([script, data, output])
            f.write(line+'\n')
    return joblist
   
if __name__ == '__main__' :
    parser = ArgumentParser()
    parser.add_argument('specie', help='name in cactus file, example: AcinonyxJubatus')
    parser.add_argument('common_name', help='name for use, example: cheetah')
    args = parser.parse_args()
    joblist = make_joblist(args.specie, args.common_name)
    #print joblist
    run_joblist(joblist, 4)

