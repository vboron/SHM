#!/usr/bin/env python3

# Using:
# agl V1.4 (c) 2020 UCL, Prof. Andrew C.R. Martin

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


def extract_data(fastadir, pdbdir):
    files = []
    error_files = []
    col = ['code' ,'VL', 'JL', 'VH', 'JH', 'angle']
    dfdata = []
    for file in os.listdir(fastadir):
        if file.endswith('.faa'):
            files.append(file)

# TODO: look at only the V sgement mutations for both (H and L)
    mismatch_data = []
    def runAGL(file, dire):
        aglresult = ''
        try:
            path = os.path.join(dire, file)
            result = (subprocess.check_output(['agl', '-a', path])).decode("utf-8")
            aglresult = result
        except subprocess.CalledProcessError:
            print(file)
            error_files.append(file)
        return aglresult

    def run_abpackingangle(pdb_code, pdb_file):
        angle = ''
        try:
            angle = (subprocess.check_output(
                ['abpackingangle', '-p', pdb_code, '-q', pdb_file])).decode("utf-8")
            angle = angle.split()
            angle = angle[1]
        except subprocess.CalledProcessError:
            error_files.append(code)
        return angle

    for file in files:
        result = runAGL(file, fastadir)
        print(result)
        result = result.replace(' ', '')
        result_data = re.split('\n', result)
        code = file[3:-4].upper()
        mismatch_data = [code]
        for data in result_data:
            if data.startswith('Mismatches'):
                data_list = data.split(':')
                mismatch = int(data_list[1])
                mismatch_data.append(mismatch)
        pdbfilepath = os.path.join(pdbdir, file[:-3]+'.pdb')
        angle = run_abpackingangle(code, pdbfilepath)
        mismatch_data.append(angle)
        dfdata.append(mismatch_data)
    df = pd.DataFrame(data=dfdata, columns=col)
    print(df)

            #     if 'Mismatches:' in data:
            #         data = data.split('Mismatches: ')
            #         mismatch_data.append([file[3:-4].upper(), data[1]])
    # df = pd.DataFrame(data=mismatch_data, columns=['code', 'mismatches'])
    # df.to_csv('agl.csv', index=False)
    # print(df)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Program for applying a rotational correction factor recursively')
    parser.add_argument(
        '--fastadir', help='Directory of fasta files', required=True)
    parser.add_argument(
        '--pdbdir', help='Directory of pdb files', required=True)
    args = parser.parse_args()

    extract_data(args.fastadir, args.pdbdir)