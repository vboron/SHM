#!/usr/bin/env python3
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


def extract_data(fastadir):
    files = os.listdir(fastadir)
    for file in files:
        num_res = run_abnum(file, fastadir)
        num_res = num_res.replace('------------------------------------------\n', '')
        num_res = num_res.strip()
        num_res_l = num_res.split('#')
        # num_res_l = num_res.split('\n')
        # del num_res_l[0]
        num_res_l = [i for i in num_res_l if i != '']
        print(num_res_l)
        for item in num_res_l:
            item_l = item.split('\n')
            del item_l[0]
        print(num_res_l)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compile.....')
    parser.add_argument(
        '--fastadir', help='Directory of fasta files', required=True)
    args = parser.parse_args()
    
    extract_data(args.fastadir)