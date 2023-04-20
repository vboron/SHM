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


def parse_agl_data(agl_out):
    agl_lines = agl_out.split('Chain type:')
    del agl_lines[0]
    agl_lines = [i.strip() for i in agl_lines]

    def process_chain(chain):
        chain_data = [l for l in agl_lines if l.startswith(chain)]
        for species in ['Homo sapiens', 'Mus musculus']:
            chain_data = [i.replace(f'{species}\n', 'splitter') for i in chain_data]
        chain_data = [i.replace('Mismatches:', 'splitter') for i in chain_data]
        chain_data = chain_data[0].split('splitter')
        chain_data = chain_data[1::2]
        chain_data = [i.strip().split('\n') for i in chain_data]
        germline = ''
        input_seq = ''
        for element in chain_data:
            germline = germline + element[2].strip()
            input_seq = input_seq + element[0].strip()
        chain_data = [list(a) for a in zip(list(input_seq), list(germline))]
        return chain_data
    
    chainl = []
    chainh = []
    if any('Light' in line for line in agl_lines):
        chainl = process_chain('Light')
    if any('Heavy' in line for line in agl_lines):
        chainh = process_chain('Heavy')
    return chainl, chainh


def label_res_mut(l_muts, h_muts, l_num, h_num):
    l_list = []
    h_list = []

    def pair_pos_num_w_res(mut_list, num_list, comb_list):
        n = 0
        m = 0
        while n < len(num_list) and m < len(mut_list):
            if num_list[n][1] == mut_list[m][0]:
                comb_list.append([num_list[n][0], mut_list[m][0], mut_list[m][1]])
                n+=1
                m+=1
            else:
                m+=1
    
    pair_pos_num_w_res(l_muts, l_num, l_list)
    pair_pos_num_w_res(h_muts, h_num, h_list)
    num_mut_pairs = l_list + h_list
    return num_mut_pairs



def find_hydrophobicity_for_positions(fastadir):
    def cal_hydrophob_change(df):
        df['input_hydrophob'] = df['input'].map(
            lambda x: utils_shm.hydrophobicity(x))
        df['germ_hydrophob'] = df['germline'].map(
            lambda x: utils_shm.hydrophobicity(x))
        change_in_hydrophobicity = df['germ_hydrophob'].sum() - df['input_hydrophob'].sum()
        return change_in_hydrophobicity
    

    def calc_hydrophobicity_for_loops(df):
        cdrL1_pos = [f'L{i}' for i in range(24, 35)]
        cdrL2_pos = [f'L{i}' for i in range(50, 57)]
        cdrL3_pos = [f'L{i}' for i in range(89, 98)]
        cdrH1_pos = [f'H{i}' for i in range(31, 36)]
        cdrH2_pos = [f'H{i}' for i in range(50, 59)]
        # cdrH3_pos = [f'H{i}' for i in range(95, 103)]

        def hydrophob_for_loop(pos_list, df):
            df_loop = df[df['L/H position'].isin(pos_list)]
            hydrophob_change = f'{cal_hydrophob_change(df_loop):.2f}'
            return hydrophob_change
        
        dH_l1 = hydrophob_for_loop(cdrL1_pos, df)
        dH_l2 = hydrophob_for_loop(cdrL2_pos, df)
        dH_l3 = hydrophob_for_loop(cdrL3_pos, df)
        dH_h1 = hydrophob_for_loop(cdrH1_pos, df)
        dH_h2 = hydrophob_for_loop(cdrH2_pos, df)
        
        return dH_l1, dH_l2, dH_l3, dH_h1, dH_h2
        
    
    files = os.listdir(fastadir)
    hydrophob_data = []
    for file in files:
        resl, resh = extract_abnum_data(file, fastadir)
        chainl_mut, chainh_mut = parse_agl_data(run_AGL(file, fastadir))
        res_pos_pairs = label_res_mut(chainl_mut, chainh_mut, resl, resh)
        posres_df = pd.DataFrame(data=res_pos_pairs, columns=['L/H position', 'input', 'germline'])
        mut_df = posres_df[posres_df['input'] != posres_df['germline']]
        delta_hydrophobicity_all = f'{cal_hydrophob_change(mut_df):.2f}'
        dh_l1, dh_l2, dh_l3, dh_h1, dh_h2 = calc_hydrophobicity_for_loops(mut_df)
        name = file[3:-4]
        data = [name, delta_hydrophobicity_all, dh_l1, dh_l2, dh_l3, dh_h1, dh_h2]
        hydrophob_data.append(data)
    df_hydroph = pd.DataFrame(data=hydrophob_data, columns=[
                      'code', 'dh_all', 'dH_L1', 'dH_L2', 'dH_L3', 'dH_H1', 'dH_H2'])
    df_hydroph = df_hydroph.astype({'dh_all': 'float64', 'dH_L1': 'float64', 'dH_L2': 'float64', 'dH_L3': 'float64', 
                                    'dH_H1': 'float64', 'dH_H2': 'float64'})
    print(df_hydroph)
    return df_hydroph


def extract_mut_data(fastadir):
    files = os.listdir(fastadir)
    col = ['code', 'VL', 'VH']
    mismatch_dic = {}
    mismatch_data = []

    for file in files:
        print(file)
        result = run_AGL(file, fastadir)
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
            if t[0].startswith('VH'):
                mismatch_s = t[-1].split(':')
                mismatch_dic['VH'] = int(mismatch_s[1])
            if t[0].startswith('VL'):
                mismatch_s = t[-1].split(':')
                mismatch_dic['VL'] = int(mismatch_s[1])
            if 'VL' not in mismatch_dic:
                mismatch_dic['VL'] = 0
            if 'VH' not in mismatch_dic:
                mismatch_dic['VH'] = 0

        print([mismatch_dic['VL'], mismatch_dic['VH']])
    #     result = result.replace(' ', '')
    #     code = file[3:-4]
    #     mismatch_data = [code]
    #     for data in result_data:
    #         if data.startswith('Mismatches'):
    #             data_list = data.split(':')
    #             mismatch = data_list[1]
    #             mismatch_data.append(mismatch)
    #     dfdata.append(mismatch_data)
    # df = pd.DataFrame(data=dfdata, columns=col)
    # df.dropna(inplace=True)
    # df = df.astype({'VH': 'int64', 'VL': 'int64'})
    # df['total_mut'] = df['VL'] + df['VH']
    return 

# TODO
# Graph dH/mutation vs count

def combine_mut_hydrophob(hydrophob_df, mut_df):
    final_df = pd.merge(mut_df, hydrophob_df, on='code')
    final_df.to_csv('hydrophobicity.csv', index=False)
    def plot_hydrophobicity(col):
        graph.hydrophobicity_vs_mutations(x_values=final_df['total_mut'], 
                                          y_values=final_df[f'dh_{col}'], 
                                          name=f'{col}_hydrophobicity')
    for value in ['all', 'l1', 'h2', 'h3']:
        plot_hydrophobicity(value)


# ********* Testing ********************************************
def run_test_parse_agl_data_bothchains():
    test_input = """
>ChainL
# Chain type: Light
VL      :  78.12% : IGKV3-20*01  : F0 : Homo sapiens
    EIVLTQ
    ||| ||
    EIVGTQ
    Mismatches: 1

JL      :  91.67% : IGKJ2*01     : F2 : Homo sapiens
    YTF
    |||
    YTF
    Mismatches: 0

>ChainH
# Chain type: Heavy
VH      :  79.59% : IGHV1-3*01   : F0 : Homo sapiens
    QVQL
    ||
    QVDK
    Mismatches: 2

JH      :  75.00% : IGHJ6*03     : F2 : Homo sapiens
    PQDN
        
    YYYY
    Mismatches: 4
    """
    chainl, chainh = parse_agl_data(test_input)
    expected_chainl = [['E', 'E'], ['I', 'I'], ['V', 'V'], ['L', 'G'], ['T', 'T'], ['Q', 'Q'], ['Y', 'Y'], ['T', 'T'], ['F', 'F']]
    expected_chainh = [['Q', 'Q'], ['V', 'V'], ['Q', 'D'], ['L', 'K'], ['P', 'Y'], ['Q', 'Y'], ['D', 'Y'], ['N', 'Y']]

    utils_shm.check_equal(chainl, expected_chainl)
    utils_shm.check_equal(chainh, expected_chainh)


def run_test_parse_agl_data_singlechainL():
    test_input = """
>ChainL
# Chain type: Light
VL      :  78.12% : IGKV3-20*01  : F0 : Homo sapiens
    EIVLTQ
    ||| ||
    EIVGTQ
    Mismatches: 1

JL      :  91.67% : IGKJ2*01     : F2 : Homo sapiens
    YTF
    |||
    YTF
    Mismatches: 0
    """
    chainl, chainh = parse_agl_data(test_input)
    expected_chainl = [['E', 'E'], ['I', 'I'], ['V', 'V'], ['L', 'G'], ['T', 'T'], ['Q', 'Q'], ['Y', 'Y'], ['T', 'T'], ['F', 'F']]
    expected_chainh = []

    utils_shm.check_equal(chainl, expected_chainl)
    utils_shm.check_equal(chainh, expected_chainh)


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
    test_l_mut = [['E', 'E'], ['I', 'I'], ['V', 'V'], ['L', 'G'], ['T', 'T'], ['Q', 'Q'], ['Y', 'Y'], ['T', 'T'], ['F', 'F']]
    test_h_num = [['H1', 'Q'], ['H2', 'V'], ['H3', 'Q'], ['H4', 'L']]
    test_h_mut = [['Q', 'Q'], ['V', 'V'], ['Q', 'D'], ['L', 'K'], ['P', 'Y'], ['Q', 'Y'], ['D', 'Y'], ['N', 'Y']]

    data = label_res_mut(test_l_mut, test_h_mut, test_l_num, test_h_num)
    expected_data = [['L1', 'E', 'E'], ['L2', 'I', 'I'], ['L3', 'V', 'V'], ['L4', 'L', 'G'], ['L5', 'T', 'T'], 
                     ['H1', 'Q', 'Q'], ['H2', 'V', 'V'], ['H3', 'Q', 'D'], ['H4', 'L', 'K']]
    utils_shm.check_equal(data, expected_data)


def run_test_label_res_mut_skippedres():
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

    run_test_parse_agl_data_bothchains()
    run_test_parse_agl_data_singlechainL()
    run_test_parse_abnum_data_bothchains()
    run_test_parse_abnum_data_singlechainH()
    run_test_label_res_mut_noresiduesskippedwithin()
    # run_test_label_res_mut_skippedres()

    df_mutations = extract_mut_data(args.fastadir)
    # df_deltahydrophobicity = find_hydrophobicity_for_positions(args.fastadir)
    # combine_mut_hydrophob(df_deltahydrophobicity, df_mutations)
