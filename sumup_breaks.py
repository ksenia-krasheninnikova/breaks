#!/usr/bin/env python

from argparse import ArgumentParser
import os
from itertools import groupby

from model import Breakpoint, parse_counted

def parse_breakpoints(path, specie):
    breaks = []
    with open(path) as f:
        for line in f:
            line = line.strip().split()
            left_t, right_t = sorted([line[8],line[10]])
            left_r, right_r = sorted([line[2],line[6]])
            breaks.append((Breakpoint(left_t,right_t,specie),Breakpoint(left_r,right_r,specie)))
    return breaks


def load_breakpoints(path, suf):
    all_breaks = {}
    for name in os.listdir(path):
        if name.split('.')[-1] == suf:
            specie = name.split('.')[0]
            breaks = parse_breakpoints(os.path.join(path,name), specie)
            all_breaks[specie] = breaks
    return all_breaks


def group_breaks(all_breaks):
    bag = []
    grouped = []
    for sp in all_breaks:
        bag += all_breaks[sp]
    keyfunc = lambda x:(x[0].left,x[0].right,x[1].left,x[1].right)
    bag = sorted(bag, key=keyfunc)
    for key, brs in groupby(bag, keyfunc):
        l = [key[0], key[1], key[2], key[3]]
        sp = []
        for b in brs:
            sp += [b[0].specie]
        sp = sorted(sp)
        grouped.append(l+[sp])
    return grouped

'''
In case there is None is coordinates (i.e.end of scaffold in some breaks)
then we can assume that this should be the same as in the similar breakpoints
'''
def merge_Nones(grouped_breaks):
    nf = filter(lambda x: 'None' in x, grouped_breaks)
    f = filter(lambda x: not 'None' in x, grouped_breaks)
    res = []
    for none in nf:
        if none[1]=='None':
           major = filter(lambda x: (none[0]==x[0] or none[0]==x[1]) and none[2]==x[2] and none[3]==x[3], f)
        elif none[0]=='None':
           major = filter(lambda x: (none[1]==x[0] or none[1]==x[1]) and none[2]==x[2] and none[3]==x[3], f)
        else:
            raise Exception('None should be in reference coords!' + str(none))
        #second check if this is a new specie
        for el in major:
            f.remove(el)
        if len(major) == 1 and len(set(none[4]+major[0][4])) == len(none)+1:
            major[0][4] += none[4]
        else:
            res += [none]
        res += major
    res += f
    return res 


def parse_counted_dir_per_species(path):
    counted = {}
    for p in os.listdir(path):
        if not 'bed' in p:
            continue
        counted[p.split('.')[0]] = parse_counted(os.path.join(path,p))
    return counted

'''
if we are able to find reference and target markers involved in a breakpoint
as single-mapped in a genome that is not included into breakpoint group
then there is an outgroup that agrees with reference and we can report the breakpoint
'''
def report_results(breaks):
    for br in breaks:
        t_left = br[0]
        t_right = br[1]
        ref_left = br[2]
        ref_right = br[3]
        for sp in counted:
            sp_vals = counted[sp]
            #skip if this is not an outgroup
            if sp in br[4]:
                continue
            if ref_left in sp_vals and ref_right in sp_vals and t_left in sp_vals and t_right in sp_vals:
                print '\t'.join(map(str, br[:4]) + ['_'.join(br[4])])
                continue


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--breaks',nargs='?',default='./breaks',help='breakpoints directory')
    parser.add_argument('--suf',nargs='?',default='breaks',help='breakpoints suffix like species.breaks')
    parser.add_argument('--hits',nargs='?',default='./hits',help='breakpoints directory with *.bed files')
    parser.add_argument('--ref')
    args = parser.parse_args()
    reference_markers = []
    breaks = load_breakpoints(args.breaks,args.suf)
    breaks = group_breaks(breaks)
    breaks = merge_Nones(breaks)
    counted = parse_counted_dir_per_species(args.hits)
    report_results(breaks)

