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
import graphing_shm as graph


def filter_line(l, free, complexed):
    l = l.strip()
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
             1mrc_0:1mrf_0H
             4org_1:"""
    free = []
    complexed = []
    for l in inp.split('\n'):
        filter_line(l, free, complexed)
    assert free == ['3u36_0,3u36_3,3u36_2,3u36_1', '1mrc_0', '4org_1']
    assert complexed == ['1zea_0PH', '1mrf_0H']
    print('test done')


def parse_redund_file(red_file):
    with open(red_file, 'r') as f:
        lines = [l.strip() for l in f.readlines()]
        free = []
        complexed = []
        all_redund = []
        for line in lines:
            filter_line(line.strip(), free, complexed)
            temp_line = line.replace(':', ',')
            all_redund.append(temp_line)
    print(all_redund)
    return free, complexed, all_redund


def dict_for_names(free, complexed, all_red):
    free_dic = {}
    complexed_dic = {}
    full_red_dic = {}
    def make_dic(list, dict):
        for item in list:
            if ',' in item:
                split_item = item.split(',')
                for si in split_item:
                    dict[si] = split_item[0]
            else:
                dict[item] = item
    make_dic(free, free_dic)
    make_dic(complexed, complexed_dic)
    make_dic(all_red, full_red_dic)
    return free_dic, complexed_dic, full_red_dic


def extract_data(fastadir, pdbdir, files, dictionary):
    error_files = []
    col = ['code', 'VL', 'JL', 'VH', 'JH', 'angle']
    dfdata = []
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
        code = file[3:-4]
        redund_code = dictionary[code]
        mismatch_data = [redund_code]
        for data in result_data:
            if data.startswith('Mismatches'):
                data_list = data.split(':')
                mismatch = int(data_list[1])
                mismatch_data.append(mismatch)
        pdbfilepath = os.path.join(pdbdir, file[:-3]+'cho')
        angle = run_abpackingangle(code.upper(), pdbfilepath)
        mismatch_data.append(angle)
        dfdata.append(mismatch_data)
    df = pd.DataFrame(data=dfdata, columns=col)
    try:
        df = df[df['angle'].str.contains('Packing') == False]
    except:
        print('No missing angles.')
    df['angle'] = df['angle'].astype(float)
    df.dropna(inplace=True)
    df.drop(columns=['JL', 'JH'])
    col.remove('JH')
    col.remove('JL')
    aggregation_func = {'angle': ['max', 'min']}
    df = df.groupby(col[:-1]).aggregate(aggregation_func)
    df['angle_range'] = df[('angle', 'max')]-df[('angle', 'min')]
    df = df.reset_index()
    df.columns = df.columns.droplevel(level=0)
    col.remove('angle')
    col = col + ['angle_min', 'angle_max', 'angle_range']
    df.columns = col
    df['total_mut'] = df['VL'] + df['VH']
    return df


def find_maxrange_per_mutation_count(df, mut_col):
    max_df = df.groupby(mut_col).angle_range.agg(['max'])
    max_df.reset_index(inplace=True)
    max_df.columns = [mut_col, 'max_angle_range']
    return max_df


def make_graphs(df, graph_name):
    cols = ['VL', 'VH', 'total_mut']
    for col in cols:
        max_df = find_maxrange_per_mutation_count(df, col)
        if col == 'total_mut':
            graph.mutations_vs_angrange(
                df, col, 'VH + VL', './', f'agl_{graph_name}_{col}_graph', max_df)
        else:
            graph.mutations_vs_angrange(
                df, col, col, './', f'agl_{graph_name}_{col}_graph', max_df)


def run_for_free_complexed(fastadir, pdbdir, free_d, complexed_d, both_d, proportion):
    files = []
    for file in os.listdir(fastadir):
        if file.endswith('.faa'):
            files.append(file)
    files.sort()

    def find_mut(files_l, dic):
        df = extract_data(fastadir, pdbdir, files_l, dic)
        df = df.sort_values(by='angle_range', ascending=False)
        top_x = len(df.index) * float(proportion)
        if top_x < 1:
            top_x = 1
        df_topx = df.head(int(top_x))
        return df, df_topx
    
    def graph_topx(dictionary, group):
        files_list = [f for f in files if f[3:-4] in dictionary]
        print(f'Finding mutations for {group} antibodies...')
        df, topx_df = find_mut(files_list, dictionary)
        df.to_csv(f'{group}_mutations.csv', index=False)
        make_graphs(topx_df, group)

    graph_topx(free_d, 'free')
    graph_topx(complexed_d, 'complex')
    graph_topx(both_d, 'complex_free')


# *************************************************************************
# *** Main program                                                      ***
# *************************************************************************
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
    parser.add_argument(
        '--top_x', help='Fraction of samples, sorted by total number of mutations, which will be graphed', required=True)
    args = parser.parse_args()

    free_list, complex_list, all_list = parse_redund_file(args.redfile)
    # dict_free, dict_complex, dict_all= dict_for_names(free_list, complex_list, all_list)
    # run_for_free_complexed(args.fastadir, args.pdbdir,
    #                        dict_free, dict_complex, dict_all, args.top_x)
