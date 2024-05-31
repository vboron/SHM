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
            ['agl', '-d', '-a', path])).decode("utf-8")
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

def extract_mut_data(fastadir):
    hydrophobicity_class = {'V': 'hydrophobic', 'I': 'hydrophobic', 'F': 'hydrophobic', 'L': 'hydrophobic', 'W': 'hydrophobic', 
                            'M': 'hydrophobic', 'R': 'hydrophilic', 'K': 'hydrophilic', 'D': 'hydrophilic', 
                            'Q': 'hydrophilic', 'N': 'hydrophilic', 'E': 'hydrophilic', 'H': 'hydrophilic', 
                            'S': 'hydrophilic'}

    def cal_hydrophob_change(df):
        df['input_class'] = df['input'].map(hydrophobicity_class)
        df['germ_class'] = df['germline'].map(hydrophobicity_class)
        df_clear = df[df['input_class'] != df['germ_class']]
        induced_hydrophilic = ''
        induced_hydrophobic = ''
        try:
            induced_hydrophilic = int(df_clear['input_class'].value_counts()['hydrophilic'])
        except:
            induced_hydrophilic = 0
        try:
            induced_hydrophobic = int(df_clear['input_class'].value_counts()['hydrophobic'])
        except:
            induced_hydrophobic = 0
        return [induced_hydrophilic, induced_hydrophobic]
    
    def calc_hydrophobicity_for_loops(df):
        print('Full df:\n', df)
        cdrL1_pos = [f'L{i}' for i in range(24, 35)]
        cdrL2_pos = [f'L{i}' for i in range(50, 57)]
        cdrL3_pos = [f'L{i}' for i in range(89, 98)]
        cdrH1_pos = [f'H{i}' for i in range(31, 36)]
        cdrH2_pos = [f'H{i}' for i in range(50, 59)]
        cdrH3_pos = [f'H{i}' for i in range(95, 103)] + [f'H100{i}' for i in char_range('A', 'K')]
        cdr_pos = cdrL1_pos + cdrL2_pos + cdrL3_pos + cdrH1_pos + cdrH2_pos + cdrH3_pos
        fwk_pos = [i for i in df['L/H position'].tolist() if i not in cdr_pos]

        def hydrophob_for_loop(pos_list):
            print(pos_list)
            df_loop = df[df['L/H position'].isin(pos_list)]
            print(df_loop)
            loop_len = len(df_loop.index)
            df_loop = df_loop[df_loop['input'] != df_loop['germline']]
            loop_mut_count = len(df_loop.index)
            
            return [cal_hydrophob_change(df_loop), loop_len, loop_mut_count]

        dH_l1_data = hydrophob_for_loop(cdrL1_pos)
        dH_l2_data = hydrophob_for_loop(cdrL2_pos)
        dH_l3_data = hydrophob_for_loop(cdrL3_pos)
        dH_h1_data = hydrophob_for_loop(cdrH1_pos)
        dH_h2_data = hydrophob_for_loop(cdrH2_pos)
        dH_h3_data = hydrophob_for_loop(cdrH3_pos)
        dH_all_loops_data = hydrophob_for_loop(cdr_pos)
        dH_fwk_data = hydrophob_for_loop(fwk_pos)
        
        return dH_l1_data, dH_l2_data, dH_l3_data, dH_h1_data, dH_h2_data, dH_h3_data, dH_all_loops_data, dH_fwk_data
    
    files = os.listdir(fastadir)
    tot_files = len(files)
    current_file = 0
    hydrophob_data = []
    abnum_fasta_dir = 'fasta_abnum'
    os.mkdir(abnum_fasta_dir)

    for file in files:
        print(file)
        resl, resh = extract_abnum_data(file, fastadir)
        abnum_seq = ''.join(['>ChainL\n'] + [l[1] for l in resl] + ['\n>ChainH\n'] + [h[1] for h in resh])
        abnum_file = f'abnum_{file}'
        abnum_file_path = os.path.join(abnum_fasta_dir, abnum_file)
        with open(abnum_file_path, 'w') as f:
            f.write(abnum_seq)

        input_L = ''
        germline_L = ''
        input_H = ''
        germline_H = ''
        result = run_AGL(abnum_file, abnum_fasta_dir)
        print(result)
        temp = result.replace('\n# ', 'splitter')
        temp = temp.replace('\n\n', 'splitter')
        temp = temp.split('splitter')
        temp = [t for t in temp if '>' not in t] 
        temp = [t.replace('Chain type: Heavy\n', '') for t in temp]
        temp = [t.replace('Chain type: Light\n', '') for t in temp]
        temp = [t for t in temp if t.startswith('VH') or t.startswith('VL')]
        temp = [t.replace(' ', '') for t in temp]
        temp = [t.split('\n') for t in temp]
        for t in temp:
            if t[0].startswith('VL'):
                input_L = input_L + t[1]
                germline_L = germline_L + t[3]
                input_L_list, germline_L_list = equalize_lists(input_L, germline_L)
            if t[0].startswith('VH'):
                input_H = input_H + t[1]
                germline_H = germline_H + t[3]
                input_H_list, germline_H_list = equalize_lists(input_H, germline_H)

        l_mut = [list(a) for a in zip(input_L_list, germline_L_list)]
        print(l_mut)
        h_mut = [list(a) for a in zip(input_H_list, germline_H_list)]
        print(h_mut)
        
        res_pos_pairs = label_res_mut(l_mut, h_mut, resl, resh)
        posres_df = pd.DataFrame(data=res_pos_pairs, columns=['L/H position', 'input', 'germline'])
        print(posres_df)
        mut_df = posres_df[posres_df['input'] != posres_df['germline']]
        mut_count = mut_df.shape[0]
        dh_all = cal_hydrophob_change(mut_df)
        data = [file[:-4], int(mut_count), dh_all[0], dh_all[1]]
        dh_l1, dh_l2, dh_l3, dh_h1, dh_h2, dh_h3, dh_cdrs, dh_fwk = calc_hydrophobicity_for_loops(posres_df)
        for region in [dh_cdrs, dh_fwk, dh_l1, dh_l2, dh_l3, dh_h1, dh_h2, dh_h3]:
            data.extend([region[0][0], region[0][1], region[1], region[2]])

        hydrophob_data.append(data)
        current_file += 1
        print(f'Progress: {current_file/tot_files*100:.2f}%')

    df_hydroph = pd.DataFrame(data=hydrophob_data, columns=[
                      'code','mut_count', 'hydrophilics_all', 'hydrophobics_all', 
                      'hydrophilics_CDRs', 'hydrophobics_CDRs', 'len_CDRs', 'mut_count_CDRs', 
                      'hydrophilics_FWk', 'hydrophobics_FWk', 'len_FWk', 'mut_count_framework', 
                      'hydrophilics_L1', 'hydrophobics_L1', 'len_L1', 'mut_count_L1',
                      'hydrophilics_L2', 'hydrophobics_L2', 'len_L2', 'mut_count_L2',
                      'hydrophilics_L3', 'hydrophobics_L3', 'len_L3', 'mut_count_L3',
                      'hydrophilics_H1', 'hydrophobics_H1', 'len_H1', 'mut_count_H1',
                      'hydrophilics_H2', 'hydrophobics_H2', 'len_H2', 'mut_count_H2',
                      'hydrophilics_H3', 'hydrophobics_H3', 'len_H3', 'mut_count_H3'])

    df_hydroph.sort_values('mut_count', inplace=True)

    df_final_hydroph = df_hydroph[2:].groupby('mut_count').aggregate('mean').reset_index()
    # print(df_final_hydroph)
    df_dist = df_hydroph[['code','mut_count', 'len_CDRs', 'hydrophilics_CDRs', 'hydrophobics_CDRs']]
    df_dist['fraction_hydrophilic'] = df_dist['hydrophilics_CDRs'] / df_dist['len_CDRs']
    df_dist['fraction_hydrophobic'] = df_dist['hydrophobics_CDRs'] / df_dist['len_CDRs']
    
    df_dist.to_csv('fractional_hydrophobicity_data.csv', index=False)
    df_final_hydroph.to_csv('introduced_hydrophobicity_data.csv', index=False)
    
    graph.introduced_hydrophobicity(df_final_hydroph)
    graph.introduced_fractional_hydrophobicity(df_dist)

    return


# ********* Testing ********************************************
def run_test_parse_abnum_data_bothchains():
    test_input = """
    # Numbered sequence  1
    L1 E
    L2 I
    L3 V
    L4 L
    L5 T
    ------------------------------------------
    # Numbered sequence  2
    H1 Q
    H2 V
    H3 Q
    H4 L
    ------------------------------------------

    """
    resnuml, resnumh = parse_abnum_data(test_input)
    expected_resnuml = [['L1', 'E'], ['L2', 'I'], ['L3', 'V'], ['L4', 'L'], ['L5', 'T']]
    expected_resnumh = [['H1', 'Q'], ['H2', 'V'], ['H3', 'Q'], ['H4', 'L']]

    utils_shm.check_equal(resnuml, expected_resnuml)
    utils_shm.check_equal(resnumh, expected_resnumh)


def run_test_parse_abnum_data_singlechainH():
    test_input = """
    # Numbered sequence  1
    H1 Q
    H2 V
    H3 Q
    H4 L
    ------------------------------------------

    """
    resnuml, resnumh = parse_abnum_data(test_input)
    expected_resnuml = []
    expected_resnumh = [['H1', 'Q'], ['H2', 'V'], ['H3', 'Q'], ['H4', 'L']]

    utils_shm.check_equal(resnuml, expected_resnuml)
    utils_shm.check_equal(resnumh, expected_resnumh)


def run_test_label_res_mut_noresiduesskippedwithin():
    test_l_num = [['L1', 'E'], ['L2', 'I'], ['L3', 'V'], ['L4', 'L'], ['L5', 'T']]
    test_l_mut = [['E', 'E'], ['I', 'I'], ['V', 'V'], ['L', 'G'], ['T', 'T']]
    test_h_num = [['H1', 'Q'], ['H2', 'V'], ['H3', 'Q'], ['H4', 'L']]
    test_h_mut = [['Q', 'Q'], ['V', 'V'], ['Q', 'D'], ['L', 'K']]

    data = label_res_mut(test_l_mut, test_h_mut, test_l_num, test_h_num)
    expected_data = [['L1', 'E', 'E'], ['L2', 'I', 'I'], ['L3', 'V', 'V'], ['L4', 'L', 'G'], ['L5', 'T', 'T'], 
                     ['H1', 'Q', 'Q'], ['H2', 'V', 'V'], ['H3', 'Q', 'D'], ['H4', 'L', 'K']]
    utils_shm.check_equal(data, expected_data)


def run_test_label_res_mut_skippedres1():
    test_l_num = []
    test_l_mut = []
    test_h_num = [['H4', 'L'], ['H97', 'P'], ['H111', 'V'], ['H112', 'S']]
    test_h_mut = [['L', 'K'], ['V', 'V'], ['S', 'Y']]

    data = label_res_mut(test_l_mut, test_h_mut, test_l_num, test_h_num)
    expected_data = [['H4', 'L', 'K'], ['H111', 'V', 'V'], ['H112', 'S', 'Y']]
    utils_shm.check_equal(data, expected_data)


def run_test_label_res_mut_skippedres2():
    test_l_num = [['L1', 'E'], ['L2', 'I'], ['L3', 'V'], ['L4', 'L'], ['L5', 'T']]
    test_l_mut = [['E', 'E'], ['I', 'I'], ['V', 'V'], ['L', 'G'], ['T', 'T'], ['Q', 'Q'], ['Y', 'Y'], ['T', 'T'], ['F', 'F']]
    test_h_num = [['H1', 'Q'], ['H2', 'V'], ['H3', 'Q'], ['H4', 'L'], ['H96', 'G'], ['H97', 'P'], ['H111', 'V'], ['H112', 'S'], ['H113', 'S']]
    test_h_mut = [['Q', 'Q'], ['V', 'V'], ['Q', 'D'], ['L', 'K'], ['V', 'V'], ['S', 'Y'], ['S', 'S']]

    data = label_res_mut(test_l_mut, test_h_mut, test_l_num, test_h_num)
    expected_data = [['L1', 'E', 'E'], ['L2', 'I', 'I'], ['L3', 'V', 'V'], ['L4', 'L', 'G'], ['L5', 'T', 'T'], 
                     ['H1', 'Q', 'Q'], ['H2', 'V', 'V'], ['H3', 'Q', 'D'], ['H4', 'L', 'K'], ['H111', 'V', 'V'], 
                     ['H112', 'S', 'Y'], ['H113', 'S', 'S']]
    utils_shm.check_equal(data, expected_data)

# *************************************************************************
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compile.....')
    parser.add_argument(
        '--fastadir', help='Directory of fasta files', required=True)
    args = parser.parse_args()

    run_test_parse_abnum_data_bothchains()
    run_test_parse_abnum_data_singlechainH()
    run_test_label_res_mut_noresiduesskippedwithin()
    # run_test_label_res_mut_skippedres1()
    # run_test_label_res_mut_skippedres2()

    extract_mut_data(args.fastadir)
