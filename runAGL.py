#!/usr/bin/env python3

# Using:
# agl V1.5 (c) 2020-22 UCL, Prof. Andrew C.R. Martin

# Usage: agl [-H|-L] [-s species] [-d datadir] [-v] [-a] [file.faa [out.txt]]
#            -H Heavy chain
#            -L Light chain
#            -s Specify a species (Homo or Mus)
#            -d Specify data directory
#            -v Verbose
#            -a Show alignments and number of mismatches

import argparse
import os
import pandas as pd
import subprocess
import re

def filter_line(l, free, complexed):
    l=l.strip()
    if '#' in l:
        print(l)
    elif l.endswith(':'):
        free.append(l[:-1])
    elif l.startswith(':'):
        complexed.append(l[1:])
    elif ':' in l:
        split_line = l.split(':')
        assert len(split_line) == 2
        free.append(split_line[0])
        complexed.append(split_line[1])

def test_filter_line():
    print('running test_filter_line...')
    inp = """:1zea_0PH
                3u36_0,3u36_3,3u36_2,3u36_1:
                1mrc_0:1mrf_0H"""
    free = []
    complexed = []
    for l in inp.split('\n'):
        filter_line(l, free, complexed)
    print(f'free={free}')
    print(f'complexed={complexed}')
    assert free == ['3u36_0,3u36_3,3u36_2,3u36_1', '1mrc_0']
    assert complexed == ['1zea_0PH', '1mrf_0H']
    print('test done')


def parse_redund_file(red_file):
    with open(red_file, 'r') as f:
        lines = [l.strip() for l in f.readlines()]
        print(len(lines))
        free = []
        complexed = []
        for line in lines:
            filter_line(line.strip(), free, complexed)
        print(len(free))
        print(len(complexed))


def extract_data(fastadir, pdbdir):
    files = []
    error_files = []
    col = ['code', 'VL', 'JL', 'VH', 'JH', 'angle']
    dfdata = []
    for file in os.listdir(fastadir):
        if file.endswith('.faa'):
            files.append(file)
    files.sort()
# TODO: look at only the V sgement mutations for both (H and L)
    mismatch_data = []

    def run_AGL(file, dire):
        aglresult = ''
        try:
            path = os.path.join(dire, file)
            result = (subprocess.check_output(
                ['agl', '-a', path])).decode("utf-8")
            aglresult = result
        except subprocess.CalledProcessError:
            print(f'AGL failed on {file}')
            error_files.append(file)
        return aglresult

    def run_abpackingangle(pdb_code, pdb_file):
        angle = 'Packing angle not found'
        try:
            angle = (subprocess.check_output(
                ['abpackingangle', '-p', pdb_code, '-q', pdb_file])).decode("utf-8")
            angle = angle.split()
            angle = angle[1]
        except subprocess.CalledProcessError:
            print(f'abpackingangle failed on {file}')
        return angle

    for file in files:
        print(file)
        result = run_AGL(file, fastadir)
        result = result.replace(' ', '')
        temp = re.split('\n', result)
        result_data = []
        for line in temp:
            if ':' in line:
                result_data.append(line)
        code = file[3:7].upper()
        mismatch_data = [code]
        for data in result_data:
            if data.startswith('Mismatches'):
                data_list = data.split(':')
                mismatch = int(data_list[1])
                mismatch_data.append(mismatch)
        pdbfilepath = os.path.join(pdbdir, file[:-3]+'cho')
        angle = run_abpackingangle(code, pdbfilepath)
        mismatch_data.append(angle)
        dfdata.append(mismatch_data)
    df = pd.DataFrame(data=dfdata, columns=col)
    try:
        df = df[df['angle'].str.contains('Packing') == False]
    except:
        print('No missing angles.')
    df.to_csv('agl_mis.csv', index=False)
    df['angle'] = df['angle'].astype(float)
    df.dropna(inplace=True)
    aggregation_func = {'angle': ['max', 'min']}
    df = df.groupby(col[:-1]).aggregate(aggregation_func)
    df['angle_range'] = df[('angle', 'max')]-df[('angle', 'min')]
    df = df.reset_index()
    df.columns = df.columns.droplevel(level=0)
    col.remove('angle')
    col = col + ['angle_min', 'angle_max', 'angle_range']
    df.columns = col
    print(df)
    df.to_csv('agl.csv', index=False)


if __name__ == '__main__':
    test_filter_line()
    parser = argparse.ArgumentParser(
        description='Compile the mutations from germline and angle ranges in redundant pdb files')
    parser.add_argument(
        '--redfile', help='List of redundant antibody files', required=True)
    parser.add_argument(
        '--fastadir', help='Directory of fasta files', required=True)
    parser.add_argument(
        '--pdbdir', help='Directory of pdb files', required=True)
    args = parser.parse_args()

    parse_redund_file(args.redfile)
    # extract_data(args.fastadir, args.pdbdir)
