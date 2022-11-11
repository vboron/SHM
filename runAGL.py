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


def extract_data(dire):
    files = []
    error_files = []
    for file in os.listdir(dire):
        if file.endswith('.faa'):
            files.append(file)


# TODO: look at only the V sgement mutations for both (H and L)
    mismatch_data = []
    for file in files:
        try:
            path = os.path.join(dire, file)
            result = (subprocess.check_output(['agl', '-a', path])).decode("utf-8")
            result_data = result.split('\n')
            for data in result_data:
                if 'Mismatches:' in data:
                    data = data.split('Mismatches: ')
                    mismatch_data.append([file[3:-4].upper(), data[1]])
        except subprocess.CalledProcessError:
            print(file)
            error_files.append(file)
    df = pd.DataFrame(data=mismatch_data, columns=['code', 'mismatches'])
    df.to_csv('agl.csv', index=False)
    print(df)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Program for applying a rotational correction factor recursively')
    parser.add_argument(
        '--dir', help='Directory with files that will be used to train models', required=True)
    args = parser.parse_args()

    extract_data(args.dir)