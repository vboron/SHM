import argparse
import os
import pandas as pd
import subprocess
import re
import graphing_shm as graph


def extract_data(fastadir):
    files = os.listdir(fastadir)
    # def run_AGL(file, dire):
    #     aglresult = ''
    #     try:
    #         path = os.path.join(dire, file)
    #         result = (subprocess.check_output(
    #             ['agl', '-a', path])).decode("utf-8")
    #         aglresult = result
    #     except subprocess.CalledProcessError:
    #         print(f'AGL failed on {file}')
    #         # error_files.append(file)
    #     return aglresult

    def run_abnum(file, dire):
        abnumresult = ''
        try:
            path = os.path.join(dire, file)
            result = (subprocess.check_output(
                ['agl', '-a', path])).decode("utf-8")
            abnumresult = result
            print(abnumresult)
        except subprocess.CalledProcessError:
            print(f'abnum failed on {file}')
            # error_files.append(file)
        return abnumresult

    for file in files:
        run_abnum(file, fastadir)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compile.....')
    parser.add_argument(
        '--fastadir', help='Directory of fasta files', required=True)
    args = parser.parse_args()
    
    extract_data(args.fastadir)