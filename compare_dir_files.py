#!/usr/bin/env python3

import argparse
import os
import shutil

def compare_dirs(dir1, dir2, new_dir_name):
    dir1_files = os.listdir(dir1)
    dir2_files = os.listdir(dir2)
    new_files = []
    for file in dir2_files:
        if file[:-4] in dir1_files:
            new_files.append(file)
    new_dir = new_dir_name
    os.mkdir(new_dir)
    print(len(new_files))
    i=0
    for file in new_files:
        i+=1
        print(i)
        print(file)
        src = os.path.join(dir2, file)
        dst = os.path.join(new_dir, file)
        shutil.copy2(src, dst)


parser = argparse.ArgumentParser(description='Program for compiling angles')
parser.add_argument('--dir1', required=True, help='Directory which is used as base')
parser.add_argument('--dir2', required=True, help='Directory from which files will be copied into new directory')
parser.add_argument('--new_dir_name', required=True, help='Name for new directory which will be created')

args = parser.parse_args()

compare_dirs(args.dir1, args.dir2, args.new_dir_name)
