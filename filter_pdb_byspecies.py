#!/usr/bin/env python3


import argparse
import os
import pandas as pd


def find_lines(dire):
    files = os.listdir(dire)
    line_terms = ['MOL_ID:', 'CHAIN:', 'ORGANISM_SCIENTIFIC:']

    for file in files:
        rel_lines = []
        with open(os.path.join(dire, file), 'r') as f:
            for line in f:
                for term in line_terms:
                    if term in line:
                        line_s = line.split(':')
                        info = line_s[1].replace(';', '')
                        rel_lines.append(info.strip())
        print(rel_lines)


# *************************************************************************
# *** Main program                                                      ***
# *************************************************************************
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compile the mutations from germline and angle ranges in redundant pdb files')
    parser.add_argument(
        '--dir', help='.', required=True)
    args = parser.parse_args()

    find_lines(args.dir)