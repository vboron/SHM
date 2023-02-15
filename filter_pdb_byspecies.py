#!/usr/bin/env python3


import argparse
import os
import shutil
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


def map_molid_chain(dire, file):
        rel_lines = []
        chain_data = {}
        with open(os.path.join(dire, file), 'r') as f:
            for line in f:
                if line.startswith('COMPND'):
                    split_line = line.split()
                    if 'MOL_ID:' in line or 'CHAIN:' in split_line[2]:
                        # print(split_line)
                        line_s = line.split(':')
                        info = line_s[1].replace(';', '')
                        rel_lines.append(info.strip())
        for x, y in pairwise(rel_lines):
            chain_data[x] = y
        return chain_data


def map_molid_org(dire, file):
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
            org_data[x] = y
        return org_data


def combine_dicts(dict1, dict2):
    dd = defaultdict(list)
    for d in (dict1, dict2): # you can list as many input dicts as you want here
        for key, value in d.items():
            dd[key].append(value)
    return dd


def map_chain_org(dire, files):
    lh_chain_species = []
    for file in files:
        print(file)
        chain_dict = map_molid_chain(dire, file)
        org_dict = map_molid_org(dire, file)
        final_dict = combine_dicts(chain_dict, org_dict)
        summary = [file[:-4]]
        l_chain = []
        h_chain = []
        for value in final_dict.values():
            if 'L' in value[0]:
                l_chain.append(value[1])
            if 'H' in value[0]:
                h_chain.append(value[1])
        if len(h_chain) == 0:
            h_chain.append(None)
        if len(l_chain) == 0:
            l_chain.append(None)
        summary = summary + h_chain + l_chain
        lh_chain_species.append(summary)
    return lh_chain_species
 

def make_df(data):
    df = pd.DataFrame(data=data, columns=['code', 'h_species', 'l_species'])
    df = df[df['h_species'].isin(['HOMO SAPIENS', 'MUS MUSCULUS'])]
    df = df[df['l_species'].isin(['HOMO SAPIENS', 'MUS MUSCULUS'])]
    human_mouse_files = df['code'].tolist()
    return human_mouse_files


def make_human_mouse_dirs(hm_list, ent_dir, pdb_dir, fasta_dir):
    def copy_files(file_type, full_dir, hm_file):
        new_dir = f'{file_type}_human_mouse'
        try:
            os.mkdir(new_dir)
        except:
            pass
        files = make_list_of_files(full_dir)
        for file in files:
            if hm_file in file:
                src = os.path.join(full_dir, file)
                dst = os.path.join(new_dir, file)
                shutil.copy2(src, dst)

    for file in hm_list:
        copy_files('ent', ent_dir, file)
        copy_files('pdb', pdb_dir, file)
        copy_files('fasta', fasta_dir, file)



# *************************************************************************
# *** Main program                                                      ***
# *************************************************************************
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compile the mutations from germline and angle ranges in redundant pdb files')
    parser.add_argument(
        '--entdir', help='.', required=True)
    parser.add_argument(
        '--pdbdir', help='.', required=True)
    parser.add_argument(
        '--fastadir', help='.', required=True)
    args = parser.parse_args()

    file_list = make_list_of_files(args.entdir)
    chain_spec = map_chain_org(args.entdir, file_list)
    hm_file_list = make_df(chain_spec)
    make_human_mouse_dirs(hm_file_list, args.entdir, args.pdbdir, args.fastadir)