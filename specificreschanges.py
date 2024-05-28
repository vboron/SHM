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
import numpy as np
from typing import List
import seaborn as sns


def run_abnum(file, dire):
    abnumresult = ''
    try:
        path = os.path.join(dire, file)
        result = (subprocess.check_output(
            ['abnum', '-f', path])).decode("utf-8")
        abnumresult = result
    except subprocess.CalledProcessError:
        print(f'abnum failed on {file}')
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
    return aglresult


def parse_abnum_data(num_res):
    split_res_num = num_res.split('\n')

    def make_LH_list(chain):
        num_res_list = [i.strip() for i in split_res_num if i.strip().startswith(chain)]
        num_res_list = [i.split(' ') for i in num_res_list]
        return num_res_list
    
    num_res_l = []
    num_res_h = []
    if any('L' in line for line in split_res_num):
        num_res_l = make_LH_list('L')
    if any('H' in line for line in split_res_num):
        num_res_h = make_LH_list('H')
    return num_res_l, num_res_h


def extract_abnum_data(in_file, dire):
    numbered_residues = run_abnum(in_file, dire)
    l_residues, h_residues = parse_abnum_data(numbered_residues)
    return l_residues, h_residues


def label_res_mut(l_muts, h_muts, l_num, h_num):
    def pair_pos_num_w_res(mut_list, num_list):
        res = []
        n = 0
        m = 0
        # if len(num_list) == len(mut_list):
        while n < len(num_list) and m < len(mut_list):
            if num_list[n][1] == mut_list[m][0]:
                # print([num_list[n][0], mut_list[m][0], mut_list[m][1]])
                res.append([num_list[n][0], mut_list[m][0], mut_list[m][1]])
                n += 1
                m += 1

        # else:
        #     'Lists are not the same lenth...'
            else:
                for m_search in range(m + 1, len(mut_list)):
                    if num_list[n][1] == mut_list[m_search][0]:
                        res.append([num_list[n][0], mut_list[m_search][0], mut_list[m_search][1]])
                        n += 1
                        m = m_search + 1
                        break
                # else executed if break never happened
                else:
                    n += 1
        return res
    
    l_list = pair_pos_num_w_res(l_muts, l_num)
    h_list = pair_pos_num_w_res(h_muts, h_num)
    return l_list + h_list

def equalize_lists(str_a, str_b, padding = 'X'):
    list_a = list(str_a)
    list_b = list(str_b)
    list_a += [padding] * (max(len(list_a), len(list_b)) - len(list_a))
    list_b += [padding] * (max(len(list_a), len(list_b)) - len(list_b))
    return list_a, list_b

def char_range(c1, c2):
    for c in range(ord(c1), ord(c2)+1):
        yield chr(c)

def dH_dy_changes(res_df):

    def cal_y_changes(df):
        # dYW_class = {'Y' : 'YW', 'W' : 'YW'}
        # df['input_class'] = df['input'].map(dYW_class)
        # df['germ_class'] = df['germline'].map(dYW_class)
        # df_clear = df[df['input_class'] != df['germ_class']]
        # germ_yw = 0
        # input_yw = 0
        # try:
        #     germ_yw = int(df_clear['germ_class'].value_counts()['YW'])
        # except:
        #     germ_yw = 0
        # try:
        #     input_yw = int(df_clear['input_class'].value_counts()['YW'])
        # except:
        #     input_yw = 0
        # return input_yw - germ_yw
        df_clear = df[df['input'] != df['germline']]
        germ_yw = 0
        input_yw = 0
        try:
            germ_yw = int(df_clear['germline'].value_counts()['Y'])
        except:
            germ_yw = 0
        try:
            input_yw = int(df_clear['input'].value_counts()['Y'])
        except:
            input_yw = 0
        return input_yw - germ_yw

    def calc_hydrophobicity_change(df):
        df['input_class'] = df['input'].map(utils_shm.hydrophobicity)
        df['germ_class'] = df['germline'].map(utils_shm.hydrophobicity)
        hydrophobiocity_change = df['input_class'].sum() - df['germ_class'].sum()
        return float(f'{hydrophobiocity_change:.2f}')

    def calc_params_for_regions(df):
        cdrL1_pos = [f'L{i}' for i in range(24, 35)]
        cdrL2_pos = [f'L{i}' for i in range(50, 57)]
        cdrL3_pos = [f'L{i}' for i in range(89, 98)]
        cdrH1_pos = [f'H{i}' for i in range(31, 36)]
        cdrH2_pos = [f'H{i}' for i in range(50, 59)]
        cdrH3_pos = [f'H{i}' for i in range(95, 103)] + [f'H100{i}' for i in char_range('A', 'K')]
        cdr_pos = cdrL1_pos + cdrL2_pos + cdrL3_pos + cdrH1_pos + cdrH2_pos + cdrH3_pos
        fwk_pos = [i for i in df['L/H position'].tolist() if i not in cdr_pos]
        fv_pos = cdr_pos + fwk_pos

        def params_for_region(pos_list):
            df_region = df[df['L/H position'].isin(pos_list)]
            region_len = len(df_region.index)
            df_region = df_region[df_region['input'] != df_region['germline']]
            loop_mut_count = len(df_region.index)
            return [cal_y_changes(df_region), calc_hydrophobicity_change(df_region), region_len, loop_mut_count]

        return params_for_region(cdr_pos) + params_for_region(fwk_pos) + params_for_region(fv_pos)
        
    return calc_params_for_regions(res_df)


def align_germline_and_get_hydrophobic_changes(fastadir):
    files = os.listdir(fastadir)
    abnum_fasta_dir = 'fasta_abnum'
    abnum_df_dir = 'abnum_csv'
    param_full_data = []

    for file in files:
        print(file)
        posres_df = pd.read_csv(os.path.join(fastadir, file))
        # resl, resh = extract_abnum_data(file, fastadir)
        # num_file_data = resl+resh
        # abnum_df = pd.DataFrame(data=num_file_data, columns=['res_pos', 'res'])
        # abnum_df.to_csv(os.path.join(abnum_df_dir, f'abnum_{file[:4]}.csv'), index=False)
        # abnum_file = f'abnum_{file}'

        # input_L = ''
        # germline_L = ''
        # input_H = ''
        # germline_H = ''
        # result = run_AGL(abnum_file, abnum_fasta_dir)
        # temp = result.replace('\n# ', 'splitter')
        # temp = temp.replace('\n\n', 'splitter')
        # temp = temp.split('splitter')
        # temp = [t for t in temp if '>' not in t] 
        # temp = [t.replace('Chain type: Heavy\n', '') for t in temp]
        # temp = [t.replace('Chain type: Light\n', '') for t in temp]
        # temp = [t for t in temp if t.startswith('VH') or t.startswith('VL')]
        # temp = [t.replace(' ', '') for t in temp]
        # temp = [t.split('\n') for t in temp]
        # for t in temp:
        #     if t[0].startswith('VL'):
        #         input_L = input_L + t[1]
        #         germline_L = germline_L + t[3]
        #         input_L_list, germline_L_list = equalize_lists(input_L, germline_L)
        #     if t[0].startswith('VH'):
        #         input_H = input_H + t[1]
        #         germline_H = germline_H + t[3]
        #         input_H_list, germline_H_list = equalize_lists(input_H, germline_H)

        # l_mut = [list(a) for a in zip(input_L_list, germline_L_list)]
        # h_mut = [list(a) for a in zip(input_H_list, germline_H_list)]
        
        # res_pos_pairs = label_res_mut(l_mut, h_mut, resl, resh)
        # posres_df = pd.DataFrame(data=res_pos_pairs, columns=['L/H position', 'input', 'germline'])
        # posres_df.to_csv(os.path.join('mutgerm_files', f'{file[:-4]}_germalign.csv'), index=False)
        param_full_data.append(dH_dy_changes(posres_df))

    cols=['cdr_dY', 'cdr_dH', 'cdr_len', 'cdr_mut', 'fwk_dY', 'fwk_dH', 'fwk_len', 'fwk_mut', 'fv_dY', 'fv_dH', 'fv_len', 'fv_mut']
    return pd.DataFrame(data=param_full_data, columns=cols)
    
def graphing_changes(param_df):
    graph.hydrophobicity_change_vs_mutation(param_df)






# *************************************************************************
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compile.....')
    parser.add_argument(
        '--fastadir', help='Directory of fasta files', required=True)
    args = parser.parse_args()

    # run_test_parse_abnum_data_bothchains()
    # run_test_parse_abnum_data_singlechainH()
    # run_test_label_res_mut_noresiduesskippedwithin()
    # run_test_label_res_mut_skippedres1()
    # run_test_label_res_mut_skippedres2()

    df_main = align_germline_and_get_hydrophobic_changes(args.fastadir)
    graphing_changes(df_main)
