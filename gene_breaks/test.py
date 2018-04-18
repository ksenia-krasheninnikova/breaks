#!/usr/bin/env python

import filecmp
import subprocess
import os
import sys


def test_all_breaks(hits, breaks, test_breaks):
    full_path_test_breaks = os.path.join(test_breaks,'case1_test.breaks')
    command = [os.path.join(get_script_path(), 'order_query_by_target.py'), '--reference', os.path.join(hits,'A.bed'), '--query', os.path.join(hits,'B.bed'),
               '| sort | uniq >', full_path_test_breaks]
    command = ' '.join(command)
    subprocess.call(command, shell=True)
    assert(filecmp.cmp(full_path_test_breaks, os.path.join(breaks,'case1.breaks'),  shallow=False)==True)


def test_seqs_breaks(hits, breaks, data, test_breaks):
    full_path_test_breaks = os.path.join(test_breaks,'case2_test.breaks')
    command = [os.path.join(get_script_path(), 'order_query_by_target.py'), '--reference', os.path.join(hits,'A.bed'), '--query', os.path.join(hits,'B.bed'),
               '--seqs', os.path.join(data, 'seqs.txt'),
               '| sort | uniq >', full_path_test_breaks]
    command = ' '.join(command)
    subprocess.call(command, shell=True)
    assert(filecmp.cmp(full_path_test_breaks, os.path.join(breaks,'case2.breaks'), shallow=False)==True)



def test_sumup_breaks(hits, breaks, data, test_breaks):
    print('Test 1...')
    test_all_breaks(hits, breaks, test_breaks)
    print('OK')
    print('Test 2...')
    test_seqs_breaks(hits, breaks, data, test_breaks)
    print('OK')
    full_path_test_sumup = os.path.join(data, 'case_sumup_test.txt')
    command = [os.path.join(get_script_path(), 'sumup_breaks.py'), '--breaks', breaks, '--hits', hits, '--ref A >', full_path_test_sumup]
    command = ' '.join(command)
    print('Test 3...')
    subprocess.call(command, shell=True)
    assert(filecmp.cmp(full_path_test_sumup, os.path.join(data,'case_sumup.txt'), shallow=False)==True)
    print('OK')


def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


if __name__ == '__main__':
    cwd = os.getcwd()
    hits = os.path.join(get_script_path(), 'data/hits/')
    breaks = os.path.join(get_script_path(), 'data/breaks')
    data = os.path.join(get_script_path(), 'data')
    test_breaks = os.path.join(cwd, 'data/test/breaks')

    os.makedirs(test_breaks)
    test_data = os.path.join(cwd, 'data')
    test_sumup_breaks(hits, breaks, data, test_breaks)
    print('Results in '+data)