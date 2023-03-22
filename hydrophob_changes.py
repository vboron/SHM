#!/usr/bin/env python3

# AGL output example (agl -a fasta_subset/pdb2ny7_0PH.faa):
#
# >ChainL
# # Chain type: Light
# VL      :  78.12% : IGKV3-20*01  : F0 : Homo sapiens
#     EIVLTQSPGTLSLSPGERATFSCRSSHSIRSRRVAWYQHKPGQAPRLVIHGVSNRASGISDRFSGSGSGTDFTLTITRVEPEDFALYYCQVYGASS
#     |||||||||||||||||||| ||| | |  |   |||| |||||||| | | | || || |||||||||||||||| | |||||| |||| || |
#     EIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSSP
#     Mismatches: 21

# JL      :  91.67% : IGKJ2*01     : F2 : Homo sapiens
#     YTFGQGTKLERK
#     |||||||||| |
#     YTFGQGTKLEIK
#     Mismatches: 1

# >ChainH
# # Chain type: Heavy
# VH      :  79.59% : IGHV1-3*01   : F0 : Homo sapiens
#     QVQLVQSGAEVKKPGASVKVSCQASGYRFSNFVIHWVRQAPGQRFEWMGWINPYNGNKEFSAKFQDRVTFTADTSANTAYMELRSLRSADTAVYYCAR
#     |||||||||||||||||||||| |||| |     |||||||||| |||||||  |||   | ||| ||| | |||| |||||| |||| |||||||||
#     QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYAMHWVRQAPGQRLEWMGWINAGNGNTKYSQKFQGRVTITRDTSASTAYMELSSLRSEDTAVYYCAR
#     Mismatches: 20

# JH      :  75.00% : IGHJ6*03     : F2 : Homo sapiens
#     PQDNYYMDVWGKGTTVIVSS
#         |||||||||||| |||
#     YYYYYYMDVWGKGTTVTVSS
#     Mismatches: 5


import argparse
import os
import pandas as pd
import subprocess
import re
import graphing_shm as graph
import utils_shm


def run_abnum(file, dire):
    abnumresult = ''
    try:
        path = os.path.join(dire, file)
        result = (subprocess.check_output(
            ['abnum', '-f', path])).decode("utf-8")
        abnumresult = result
        # print(abnumresult)
    except subprocess.CalledProcessError:
        print(f'abnum failed on {file}')
        # error_files.append(file)
    return abnumresult


def run_AGL(file, dire):
    aglresult = ''
    try:
        path = os.path.join(dire, file)
        result = (subprocess.check_output(
            ['agl', '-a', path])).decode("utf-8")
        aglresult = result
    except subprocess.CalledProcessError:
        print(f'AGL failed on {file}')
        # error_files.append(file)
    return aglresult


def parse_abnum_data(in_file, dire):
    num_res = run_abnum(in_file, dire)
    num_res = num_res.replace(
        '------------------------------------------\n', '')
    num_res = num_res.strip()
    num_res_l = num_res.split('#')
    num_res_l = [i for i in num_res_l if i != '']
    num_res_l_clean = []
    for item in num_res_l:
        item_l = item.split('\n')
        del item_l[0]
        num_res_l_clean.append(item_l)
    # print(f'{file}: {len(num_res_l_clean)}')
    all_chain_single_list = []
    for num_chain in num_res_l_clean:
        all_chain_single_list = all_chain_single_list + num_chain
    final_num_res = []
    for res in all_chain_single_list:
        res_split = res.split(' ')
        final_num_res.append(res_split)
    return final_num_res


def parse_agl_data(agl_out):
    if 'Homo sapiens' in agl_out:
        agl_out = agl_out.replace('Homo sapiens\n', 'splitter')
    if 'Mus musculus' in agl_out:
        agl_out = agl_out.replace('Mus musculus\n', 'splitter')
    agl_out = agl_out.replace('Mismatches:', 'splitter')
    agl_list = agl_out.split('splitter')
    agl_list = agl_list[1::2]
    agl_list_stripped = []
    for i in agl_list:
        i = i.strip()
        i_clean = []
        i_split = i.split('\n')
        i_split = i_split[0::2]
        for ii in i_split:
            ii = ii.strip()
            i_clean.append(ii)
        agl_list_stripped.append(i_clean)
    input_seq = ''
    germline_seq = ''
    for seq in agl_list_stripped:
        input_seq = input_seq + seq[0]
        germline_seq = germline_seq + seq[1]
    input_res_list = list(input_seq)
    germline_res_list = list(input_seq)
    input_germ_list = [list(a) for a in zip(input_res_list, germline_res_list)]
    return input_germ_list


def cal_hydrophob_change(imput_germ_pairs):
    df = pd.DataFrame(data=imput_germ_pairs, columns=['input', 'germline'])
    print(df)
    df = df[df['input'] != df['germline']]
    print(df)

def extract_data(fastadir):
    files = os.listdir(fastadir)
    for file in files:
        # numbered_res_list = parse_abnum_data(file, fastadir)
        agl_output = run_AGL(file, fastadir)
        in_germ_res_pairs = parse_agl_data(agl_output)
        print(file)
        # print(len(numbered_res_list))
        print(in_germ_res_pairs)
        cal_hydrophob_change(in_germ_res_pairs)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compile.....')
    parser.add_argument(
        '--fastadir', help='Directory of fasta files', required=True)
    args = parser.parse_args()

    extract_data(args.fastadir)
