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
        num_res_list = [i for i in split_res_num if i.startswith(chain)]
        # num_res_list = [i for i in split_res_num if i.startswith('L') or i.startswith('H')]
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
    return parse_abnum_data(numbered_residues)


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
        # print(f'Chain L: {chainl}')
    if any('Heavy' in line for line in agl_lines):
        chainh = process_chain('Heavy')
        # print(f'Chain H: {chainh}')
    return chainl, chainh


def cal_hydrophob_change(imput_germ_pairs):
    df = pd.DataFrame(data=imput_germ_pairs, columns=['input', 'germline'])
    df = df[df['input'] != df['germline']]
    df['input_hydrophob'] = df['input'].map(
        lambda x: utils_shm.hydrophobicity(x))
    df['germ_hydrophob'] = df['germline'].map(
        lambda x: utils_shm.hydrophobicity(x))
    change_in_hydrophobicity = df['germ_hydrophob'].sum() - df['input_hydrophob'].sum()
    return change_in_hydrophobicity


def extract_hydrophob_data(fastadir):
    files = os.listdir(fastadir)
    hydrophob_data = []
    for file in files:
        agl_output = run_AGL(file, fastadir)
        in_germ_res_pairs = parse_agl_data(agl_output)
        delta_hydrophobicity = f'{cal_hydrophob_change(in_germ_res_pairs):.2f}'
        name = file[3:-4]
        data = [name, delta_hydrophobicity]
        hydrophob_data.append(data)
    df = pd.DataFrame(data=hydrophob_data, columns=[
                      'code', 'delta_hydrophobicity'])
    df = df.astype({'delta_hydrophobicity': 'float64'})
    return df


def extract_mut_data(fastadir):
    files = os.listdir(fastadir)
    col = ['code', 'VL', 'JL', 'VH', 'JH']
    dfdata = []
    mismatch_data = []

    for file in files:
        print(file)
        # result = run_AGL(file, fastadir)
        # result = result.replace(' ', '')
        # temp = re.split('\n', result)
        # result_data = [l for l in temp if ':' in l]
        # code = file[3:-4]
        # mismatch_data = [code]
        # for data in result_data:
        #     if data.startswith('Mismatches'):
        #         data_list = data.split(':')
        #         mismatch = data_list[1]
        #         mismatch_data.append(mismatch)
        # dfdata.append(mismatch_data)
        print(f'abnum: {extract_abnum_data(file, fastadir)}')
        parse_agl_data(run_AGL(file, fastadir))
        # print(f'len agl: {len(parse_agl_data(run_AGL(file, fastadir)))}')
        # print([list(a) for a in zip(extract_abnum_data(file, fastadir), parse_agl_data(run_AGL(file, fastadir)))])
    # df = pd.DataFrame(data=dfdata, columns=col)
    # df.dropna(inplace=True)
    # df.drop(columns=['JL', 'JH'], inplace=True)
    # df = df.astype({'VH': 'int64', 'VL': 'int64'})
    # df['total_mut'] = df['VL'] + df['VH']
    # return df
    return


def combine_mut_hydrophob(hydrophob_df, mut_df):
    final_df = pd.merge(mut_df, hydrophob_df, on='code')
    graph.hydrophobicity_vs_mutations(
        x_values=final_df['total_mut'], y_values=final_df['delta_hydrophobicity'], name='hydrophobicity')


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

    # df_deltahydrophobicity = extract_hydrophob_data(args.fastadir)
    df_mutations = extract_mut_data(args.fastadir)
    # combine_mut_hydrophob(df_deltahydrophobicity, df_mutations)

#####TODO
#Remove the blank [''] from the abnum list