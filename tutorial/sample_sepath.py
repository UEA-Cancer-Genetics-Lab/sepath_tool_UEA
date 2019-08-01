#!/usr/bin/env python

#This is a sample script that will download some useful files and launch a test version of SEPATH

#import some basic modules
import os
import glob
import re
import sys
import subprocess

#print welcome message
print('\033[1;32;40m This is a sample script that will download some useful files and launch a test version of SEPATH\n\033[0;37;40m')

#run a dependencies check
print('\033[1;32;40m\n---Checking dependencies---\n\033[0;37;40m')
#set list of tools required
error_tools = []

try:
    import pysam
except:
    print('\033[1;31;40m Pysam not found - Python module required if working from BAM input\033[0;37;40m')
    error_found=True

#list of system commands to try
commands = [
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
        error_tools.append('\033[1;31;40m %s not found - Please install and check functionality\033[0;37;40m' %tool)
        error_found=True

if (len(error_tools) > 0):
    print('\033[1;31;40m Errors found: \033[0;37;40m' )
    for tool in error_tools:
        print(tool)
    print('\nCheck installation instructions for help installing each tool')
    raise ValueError('\033[1;31;40m -- Dependencies not in working order... Exiting Script \033[0;37;40m')
else:
    print('\033[1;32;40m All software dependencies seem in order!\n\033[0;37;40m')


# Check if dmp/ directory exists - if not prompt user to download some useful files
if (not os.path.isdir('dmp')) or (not os.path.exists('dmp/minikraken_20171101_4GB_dustmasked/database.kdb')) or (not os.path.exists('dmp/grch38.fna')):
    print('\033[1;32;40m\n---Creating dmp directory---\033[0;37;40m \n')
    os.system('mkdir -pv dmp')
    print('\033[1;32;40m\n---Downloading Mini-Kraken Database--- \033[0;37;40m \n')
    os.system('wget -c https://ccb.jhu.edu/software/kraken/dl/minikraken_20171101_4GB_dustmasked.tgz -O dmp/minikraken_dust.tgz')
    os.system('tar -xzvf dmp/minikraken_dust.tgz; mv minikraken_20171101_4GB_dustmasked/ dmp/')
    print('\033[1;32;40m \n---Downloading Sample human reference---\033[0;37;40m \n')
    os.system('wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.fna.gz -O dmp/grch38.fna.gz; gunzip dmp/grch38.fna.gz; head -n 4000 dmp/grch38.fna > dmp/small_human.fna')
    # Check files downloaded successfully
    if (not os.path.isdir('dmp')) or (not os.path.exists('dmp/minikraken_20171101_4GB_dustmasked/database.kdb')) or (not os.path.exists('dmp/grch38.fna')):
        raise ValueError('\033[1;31;40m \nSomething went wrong during the download of minikraken or human reference databases - Please debug and re-run script and submit a bug report \033[0;37;40m \n')
else:
    print('\033[1;32;40m\n---Mini-Kraken and Human Reference Files Located Successfully---\033[0;37;40m \n')

#Check that bam/BAM files exist in the files/ directory
bam_files = glob.glob('files/*.BAM')
if len(bam_files) == 0:
    bam_files = glob.glob('files/*.bam')

#check that bam files exist
if len(bam_files) == 0:
    raise ValueError('\033[1;31;40m \n No BAMfiles found in the files/ directory -- please put some in there or re-download the tutorial ones \033[0;37;40m \n')
else:
    print('\033[1;32;40m \n%s BAM files found in the files/ directory\033[0;37;40m \n' %len(bam_files))


# Launch a complete snakemake and make directories required
print('\033[1;32;40m\n---Launching Snakemake---\033[0;37;40m \n')
os.system('mkdir -pv output')
os.system('snakemake -s sepath.snake --keep-going')

