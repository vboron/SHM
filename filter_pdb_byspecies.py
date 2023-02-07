#!/usr/bin/env python3


import argparse
import os
import pandas as pd


def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)


def find_lines(dire):
    files = os.listdir(dire)
    line_terms = ['MOL_ID:', 'CHAIN:', 'ORGANISM_SCIENTIFIC:']
    species_info = []
    for file in files:
        l_chain = []
        h_chain = []
        summary = [file[:-4]]
        data = {}
        rel_lines = []
        with open(os.path.join(dire, file), 'r') as f:
            for line in f:
                for term in line_terms:
                    if term in line:
                        line_s = line.split(':')
                        info = line_s[1].replace(';', '')
                        rel_lines.append(info.strip())
        for x, y in pairwise(rel_lines):
            if x in data:
               data[x].append(y) 
            else:
                data[x] = [y]
        for value in data.values():
            if 'L' in value[0]:
                # l_chain.append('L')
                l_chain.append(value[1])
            if 'H' in value[0]:
                # h_chain.append('H')
                h_chain.append(value[1])
        summary = summary + h_chain + l_chain
        species_info.append(summary)
    return species_info


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

    chain_spec = find_lines(args.dir)
    make_df(chain_spec)