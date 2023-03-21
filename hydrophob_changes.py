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


def extract_data(fastadir):
    files = os.listdir(fastadir)
    for file in files:
        # num_res = run_abnum(file, fastadir)
        # num_res = num_res.replace('------------------------------------------\n', '')
        # num_res = num_res.strip()
        # num_res_l = num_res.split('#')
        # num_res_l = [i for i in num_res_l if i != '']
        # num_res_l_clean = []
        # for item in num_res_l:
        #     item_l = item.split('\n')
        #     del item_l[0]
        #     num_res_l_clean.append(item_l)
        # print(f'{file}: {len(num_res_l_clean)}')

        agl_out = run_AGL(file, fastadir)
        if 'Homo sapiens' in agl_out:
            agl_out = agl_out.replace('Homo sapiens\n', 'h_or_m')
        if 'Mus musculus' in agl_out:
            agl_out = agl_out.replace('Mus musculus\n', 'h_or_m')
    
        agl_list = agl_out.split('h_or_m')
        # del agl_list[0]
        # print(agl_list)
        agl_list_split = []
        for i in agl_list:
            i = i.strip()
            i = i.split('Mismatches:')
            # i = i.strip()
            agl_list_split.append(i)
        print(file)
        print(agl_list_split)
        # print(f'{file}: {len(agl_list)}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compile.....')
    parser.add_argument(
        '--fastadir', help='Directory of fasta files', required=True)
    args = parser.parse_args()
    
    extract_data(args.fastadir)