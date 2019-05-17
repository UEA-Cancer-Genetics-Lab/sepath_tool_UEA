#!/usr/bin/env python3

# This script is used to check dependencies for SEPATH

import os
import subprocess
import sys

print('Checking SEPATH dependencies\n')

#set list of tools required
error_tools = []

try:
    import pysam
except:
    print('Pysam not found - Python module required if working from BAM input')

#list of system commands to try
commands = [
    'SEPATH.py --version'
    'snakemake --version',
    'java -version',
    'trimmomatic -version',
    'bbduk.sh -version',
    'motus --version',
    'kraken -version',
    'metaspades.py --version'
]

#try each command in the shell - suppress output by sending to devnull
for command in commands:
    tool = command.split(' ')[0]
    out_null = open(os.devnull, 'w')
    try:
        subprocess.check_call(command, stdout=out_null, stderr=out_null, shell=True)
    except:
        error_tools.append('%s not found - Please install and check functionality' %tool)

#Set some warning colours to highlight dependencies missing
class warning_cols:
    WARNING = '\x1b[1;31m'
    RESET = '\033[0m'

if (len(error_tools) > 0):
    print('%sErrors found:%s' %(warning_cols.WARNING, warning_cols.RESET))
    for tool in error_tools:
        print(tool)
    print('\nCheck installation instructions for help installing each tool')
else:
    print('All seems to be in working order!')
