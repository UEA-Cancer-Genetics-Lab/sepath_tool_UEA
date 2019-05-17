#!/usr/bin/env python

#This script will launch snakemake commands for paired end fastq filtering and assembly

import os
import glob
import time
import argparse

#build argument parser
parser = argparse.ArgumentParser(description='Searching for Pathogens - SEPATH - script to run quality trimming/host depletion and assembly on paired end fastq input')
parser.add_argument('--outdir', nargs=1, required=True, help='enter the absolute path of your output directory')

#parse arguments and set as variables
args = parser.parse_args()
cwd = os.getcwd()
infiles = ('%s/input_files.txt' %cwd)
outdir = str(args.outdir[0])

#load up that input file list:
with open(infiles) as f:
	input_files = f.read().splitlines()

#sort to make sure they are in order
input_files.sort()

x=0
while x<len(input_files):
	bamfile = str(input_files[x])
	launchsnake = ('bsub -J SEPATH -q long-eth -M 10000 -R "rusage[mem=10000]" "snakemake -s  all_SEPATH.snake  --config outdir=%s input_bam=%s --cluster \\"bsub -q {params.queue} {params.cluster}\\" --jobs 20 --latency-wait 1200 --timestamp"' %(outdir, bamfile))
	#develop ID to look for from snakemake output - or just keep it as jobdone for now?
	FTL = 'jobdone'
	#launch the snakemake
	os.system(launchsnake)
	while(os.path.isfile(FTL) == False):
		print('Job not finished, Checking again in 10 minutes')
		time.sleep(600)
	#do some house keeping
	os.system('rm -rf *fq.gz jobdone')
	x=x+1







