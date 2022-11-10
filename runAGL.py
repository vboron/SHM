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

    for file in files:
        try:
            path = os.path.join(dire, file)
            result = (subprocess.check_output(['agl', '-H|-L', '-v', '-a', path])).decode("utf-8")
            print(result)
        except subprocess.CalledProcessError:
            error_files.append(file)
# angle = (subprocess.check_output(
#                 ['abpackingangle', '-p', pdb_code, '-q', pdb_file])).decode("utf-8")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Program for applying a rotational correction factor recursively')
    parser.add_argument(
        '--dir', help='Directory with files that will be used to train models', required=True)
    args = parser.parse_args()

    extract_data(args.dir)