#!/usr/bin/env python3
import os
import sys
import shutil
#
cwd = os.getcwd()
for file in os.listdir(sys.argv[1]):
    if file.endswith('.mar'):
        src = os.path.join(cwd, sys.argv[1], file)
        dst = os.path.join(cwd, sys.argv[1], (file[3:9].upper() + '.pdb'))
        # os.rename(src, dst)
        shutil.copy2(src, dst)
    # if file.endswith('.pdb'):
    #     src = os.path.join(cwd, sys.argv[1], file)
    #     dst = os.path.join(cwd, sys.argv[1], (file[1:6] + '.pdb'))
    #     os.rename(src, dst)
