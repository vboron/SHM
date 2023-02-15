#!/usr/bin/env python3


import argparse
import os
import pandas as pd
from collections import defaultdict


def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)

def combine_dict(d1, d2):
    return {k: tuple(d[k] for d in (d1, d2) if k in d) for k in set(d1.keys()) | set(d2.keys())}

def make_list_of_files(dire):
    files = os.listdir(dire)
    return files

def map_molid_chain(dire, files):
    for file in files:
        rel_lines = []
        summary = [file[:-4]]
        l_chain = []
        h_chain = []
        chain_data = {}
        with open(os.path.join(dire, file), 'r') as f:
            for line in f:
                if line.startswith('COMPND'):
                    if 'MOL_ID:' in line or 'CHAIN:' in line:
                        line_s = line.split(':')
                        info = line_s[1].replace(';', '')
                        rel_lines.append(info.strip())
        # print(rel_lines)
        for x, y in pairwise(rel_lines):
            if 'H' in y or 'L' in y:
                if x in chain_data:
                    chain_data[x].append(y) 
                else:
                    chain_data[x] = [y]
        prev = ''
        org_list = []
        with open(os.path.join(dire, file), 'r') as f:
            for line in f:
                if line.startswith('SOURCE'):
                    for term in ['MOL_ID:', 'ORGANISM_SCIENTIFIC:']:
                        if term in line:
                            # print(line)
                            if prev == term:
                                # print('missing line...')
                                org_list.append(None)
                            line_s = line.split(':')
                            info = line_s[1].replace(';', '').strip()
                            org_list.append(info)
                            prev = term
            if prev == 'MOL_ID:':
                org_list.append(None)
        # print(org_list)
        org_data = {}
        for x, y in pairwise(org_list):
            if x in org_data:
                org_data[x].append(y) 
            else:
                org_data[x] = [y]
        # print(org_data)
        dd = defaultdict(list)
        for d in (chain_data, org_data): # you can list as many input dicts as you want here
            for key, value in d.items():
                dd[key].append(value)
        print(dd)
        # for value in chain_data.values():
        #     print(file, value)
        #     if 'L' in value[0]:
        #         # l_chain.append('L')
        #         l_chain.append(value[1])
        #     if 'H' in value[0]:
        #         # h_chain.append('H')
        #         h_chain.append(value[1])
        # summary = summary + h_chain + l_chain
        # print(summary)
 


def make_df(data):
    df = pd.DataFrame(data=data, columns=['code', 'h_species', 'l_species'])
    print(df)


# *************************************************************************
# *** Main program                                                      ***
# *************************************************************************
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compile the mutations from germline and angle ranges in redundant pdb files')
    parser.add_argument(
        '--dir', help='.', required=True)
    args = parser.parse_args()

    file_list = make_list_of_files(args.dir)
    chain_spec = map_molid_chain(args.dir, file_list)
    # make_df(chain_spec)