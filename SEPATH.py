#!/usr/bin/env python3

import os
import sys
import argparse

# The main SEPATH overlord script - use to launch bam filtering, paried fastq filtering, bacterial classificaition or assembly/classification

# Build an argument parser
parser = argparse.ArgumentParser(description='Searching for Pathogens - SEPATH')
parser.add_argument('--mode', nargs=1, required=True, choices=['bam_filter','fastq_bacterial', 'fastq_viral', 'bam_bacterial', 'bam_viral', 'bam_all'], help='Selects the mode used to run SEPATH (bam_filter, fastq_filter, bam_bacterial, bam_viral, bam_all)')
parser.add_argument('--num_dirs', nargs=1, type=int, help='Enter number of working directories', default=1)
parser.add_argument('--in_files', nargs=1, help='Enter list of inputfiles (default input_files.txt)', default='input_files.txt')
#parser.add_argument('--files_per_sample', nargs=1, type=int, required=True, help='Enter the number of files per sample, ie BAM files = 1, paired end = 2, filtered paired end = 3 (one single end)')
parser.add_argument('--outdir', nargs=1, required=True, help='enter the absolute path of your output directory')

#version
parser.add_argument('--version', action='version', version='SEPATH Version 1.0 BETA', help='Returns the version of SEPATH')

args = parser.parse_args()

#Declare SEPATH run mode
if (args.mode[0] == 'bam_filter'):
    print('Running BAM filtering pipeline')
elif (args.mode[0] == 'fastq_bacterial'):
    print('Running paired end fastq filter pipeline')
elif (args.mode[0] == 'fastq_viral'):
    print('Running Fastq metagenomic assembly')
elif (args.mode[0] == 'bam_bacterial'):
    print('Running mOTUs2 bacterial classification pipeline')
elif (args.mode[0] == 'bam_viral'):
    print('Running metagenomic assembly & Kraken classification pipeline')
elif (args.mode[0] == 'bam_all'):
    print('Running Kraken & mOTUs2 on unmapped reads and assembled contigs')

#set command to launch jobs for your cluster
long_job = 'bsub -q long-ib -M 10000 -R "rusage[mem=10000]" -J long '

#obtain current working directory
cwd = os.getcwd()

#get location of SEPATH script to copy scripts from bin
script_path = (os.path.realpath(__file__)).split('SEPATH.py')[0]

#check requested output directory exists, if not, create it
outdir = str(args.outdir[0])
print("Output Directory: %s" %outdir)

if os.path.isdir(outdir) == True:
    print('Output directory: %s' %outdir)
else:
    print('Output directory not detected, creating directory: %s' %outdir)
    command = ('mkdir -pv %s' %outdir)
    os.system(command)

#declare input file list
if (args.in_files == 'input_files.txt'):
    infile = ('%s/input_files.txt' %cwd)
    print('Input file list: input_files.txt')
else:
    infile=('%s/%s' %(cwd,args.in_files[0]))
    print('Input file list: %s' %infile)

#open input file list ready for distribution into working directories
with open(str(infile)) as f:
    inputfiles = f.read().splitlines()

#ensure files are sorted by name
inputfiles.sort()
print('%s files in input file list' %len(inputfiles))
#ensure no duplicate files in input file list
for file in inputfiles:
    if inputfiles.count(file) > 1:
        raise ValueError('Duplicate files found in input file list: %s' %file)


#create working directories
work_dir_list = []
full_dir_list = []

for dir in range(1,(args.num_dirs[0]+1)):
    directory = ('workdir%s' %dir)
    full_directory = ('%s/%s' %(cwd, directory))
    command = ('mkdir -v %s' %full_directory)
    # check directories don't already exist, if it does, ruthlessly delete them and remake them
    if os.path.isdir(full_directory) == True:
        folders = ('%s already exists, would you like to overwrite it? Y/N: ' %full_directory)
        obliterate_folders = input(folders)
        if obliterate_folders.upper() == 'Y':
            os.system('rm -rf %s' %full_directory)
        else:
            raise ValueError('Overwrite permissions not provided, sort out your working directories')
    #make directory structure
    os.system(command)
    #copy all dependencies to working directories
    copy_command = ('cp %sbin/*/* %s/'   %(script_path, full_directory))
    os.system(copy_command)
    work_dir_list.append(directory)
    full_dir_list.append(full_directory)


#define a function to disperse files throughout working directories
def distribute_files(input_file_list, working_dirs, files_per_dir):
    #make sure number of input files is divisible by files per directory
    remainder = int(len(input_file_list) % files_per_dir)
    if (remainder !=0):
        raise ValueError('Number of input files not divisible by number of files per directory - you have some missing files, please check your input file list before running SEPATH again')
    #if less input files than number of directories, just print one file per directory for the number of files
    if (len(input_file_list) < len(full_dir_list)):
        counter = 0
        while (counter < len(input_file_list)):
            command = ('echo %s > %s/input_files.txt' %(input_file_list[counter], full_dir_list[counter]))
            os.system(command)
            counter = counter + 1
    #if only one directory - disperse all input files in one directory
    elif (len(full_dir_list) == 1):
        for file in input_file_list:
            command = ('echo %s >> %s/input_files.txt' %(file, full_dir_list[0]))
            os.system(command)
    else:
        #Loop through directories assigning files
        file_counter = 0
        directory_counter = 0
        while (file_counter < len(input_file_list)):
            command = ('echo %s >> %s/input_files.txt' %(input_file_list[file_counter], full_dir_list[directory_counter]))
            os.system(command)
            file_counter = file_counter + 1
            if (files_per_dir > 1):
                command = ('echo %s >> %s/input_files.txt' %(input_file_list[file_counter], full_dir_list[directory_counter]))
                os.system(command)
                file_counter = file_counter + 1
            if (files_per_dir > 2):
                command = ('echo %s >> %s/input_files.txt' %(input_file_list[file_counter], full_dir_list[directory_counter]))
                os.system(command)
                file_counter = file_counter + 1
            #reset directory counter if it has reached the end of the list of folders to ensure even distribution
            if (directory_counter == (len(full_dir_list)-1)):
                directory_counter = 0
            else:
                directory_counter = directory_counter + 1

#define function to check input files exist:
def do_files_exist(input_file_list):
    for file in input_file_list:
        if (os.path.isfile(file) == False):
            raise ValueError('File in input file list not found:\n%s' %file)
#check files exist - also written as a function in case needed later
do_files_exist(inputfiles)

#run bam filtering:
if (args.mode[0] == 'bam_filter'):
    print('Entering Bam Filtering Mode')
    #check 1 file per sample - because they're bams
    #if (int(args.files_per_sample[0]) != 1):
    #    raise ValueError('%s files per sample specified, this is not possible for BAM file processing' %(args.files_per_sample))
    #check all files end in .bam
    for file in inputfiles:
        if (file[-4:].lower() != '.bam'):
            raise ValueError('File detected in input file list that is not a bam file, please clean the input file list:\n%s' %file) 
    #distribute files among working directories
    distribute_files(inputfiles, full_dir_list, 1)# args.files_per_sample[0] changed to 1 for bam files
    #navigate to each working directory and unleash ouroboros (launch seperate snakemake job launchers) 
    # in the case of bam filtering, it is not snakemake based, so just launch the job using the long job
    for directory in full_dir_list:
        os.chdir(directory)
        #make scripts executable
        os.system('chmod +x bam_ouroboros.py')
        #develop command to unleash ouroboros upon the world (launch bam filtering job)
        bam_command = ('%s "./bam_ouroboros.py --outdir %s"' %(long_job, outdir))
        os.system(bam_command)
    #revert back to original working directory
    os.chdir(cwd)


#run fastq filtering
if (args.mode[0] == 'fastq_bacterial'):
    print('Entering FASTQ --> unmapped reads --> mOTUs2 Bacterial mode')
    #check 2 files per sample because they are paried end fastq files
    #if (int(args.files_per_sample[0]) != 2):
    #    raise ValueError('%s files per sample specified, this is not possible for paired end FASTQ processing' %(args.files_per_sample))
    #check name format is sample.fq.gz
    for file in inputfiles:
        if (file[-6:].lower() != '.fq.gz'):
            raise ValueError('Files must be in sample1_R1.fq.gz, sample1_R2.fq.gz format')
    #distribute files among working directories
    distribute_files(inputfiles, full_dir_list, 1)# args.files_per_sample[0] changed to 1 for bam files
    #navigate to each working directory and unleash ouroboros (launch seperate snakemake job launchers)
    for directory in full_dir_list:
        os.chdir(directory)
        #make scripts executable
        os.system('chmod +x fastq_ouroboros.py')
        #develop command to unleash ouroborus upon this world (launch fastq filtering job)
        fastq_command = ('%s "./fastq_ouroboros.py --outdir %s"' %(long_job, outdir))
        os.system(fastq_command)
    #revert back to original working directory
    os.chdir(cwd)


if (args.mode[0] == 'fastq_viral'):
    print('Entering FASTQ --> contig assembly --> kraken mode')
    #check 2 files per sample because they are paried end fastq files
    #if (int(args.files_per_sample[0]) != 2):
    #    raise ValueError('%s files per sample specified, this is not possible for paired end FASTQ processing' %(args.files_per_sample))
    #check name format is sample.fq.gz
    for file in inputfiles:
        if (file[-6:].lower() != '.fq.gz'):
            raise ValueError('Files must be in sample1_R1.fq.gz, sample1_R2.fq.gz format')
    #distribute files among working directories
    distribute_files(inputfiles, full_dir_list, args.files_per_sample[0])
    #navigate to each working directory and unleash ouroboros (launch seperate snakemake job launchers)
    for directory in full_dir_list:
        os.chdir(directory)
        #make scripts executable
        os.system('chmod +x fastq_assembly_ouroboros.py')
        #develop command to unleash ouroborus upon this world (launch fastq filtering job)
        fastq_command = ('%s "./fastq_assembly_ouroboros.py --outdir %s"' %(long_job, outdir))
        os.system(fastq_command)
    #revert back to original working directory
    os.chdir(cwd)


#bam --> bacteiral with retaining unmapped reads
if (args.mode[0] == 'bam_bacterial'):
    print('Entering BAM  --> Bacterial classification mode')
    #check 2 files per sample because they are paried end fastq files
    #if (int(args.files_per_sample[0]) != 1):
    #    raise ValueError('%s files per sample specified, this is not possible for paired end FASTQ processing' %(args.files_per_sample))
    #check name format is sample.fq.gz
    for file in inputfiles:
        if (file[-4:].lower() != '.bam'):
            raise ValueError('Files must be in sample1.bam format')
    #distribute files among working directories
    distribute_files(inputfiles, full_dir_list, 1)# args.files_per_sample[0] changed to 1 for bam files
    #navigate to each working directory and unleash ouroboros (launch seperate snakemake job launchers)
    for directory in full_dir_list:
        os.chdir(directory)
        #make scripts executable
        os.system('chmod +x bacterial_ouroboros.py')
        #develop command to unleash ouroborus upon this world (launch fastq filtering job)
        fastq_command = ('%s "./bacterial_ouroboros.py --outdir %s"' %(long_job, outdir))
        os.system(fastq_command)
    #revert back to original working directory
    os.chdir(cwd)


#bam --> viral with assembly mode
if (args.mode[0] == 'bam_viral'):
    print('Entering BAM  --> viral classification mode')
    #check 2 files per sample because they are paried end fastq files
    #if (int(args.files_per_sample[0]) != 1):
    #    raise ValueError('%s files per sample specified, this is not possible for paired end FASTQ processing' %(args.files_per_sample))
    #check name format is sample.fq.gz
    for file in inputfiles:
        if (file[-4:].lower() != '.bam'):
            raise ValueError('Files must be in sample1.bam format')
    #distribute files among working directories
    distribute_files(inputfiles, full_dir_list, 1)# args.files_per_sample[0] changed to 1 for bam files
    #navigate to each working directory and unleash ouroboros (launch seperate snakemake job launchers)
    for directory in full_dir_list:
        os.chdir(directory)
        #make scripts executable
        os.system('chmod +x viral_ouroboros.py')
        #develop command to unleash ouroborus upon this world (launch fastq filtering job)
        fastq_command = ('%s "./viral_ouroboros.py --outdir %s"' %(long_job, outdir))
        os.system(fastq_command)
    #revert back to original working directory
    os.chdir(cwd)

#bam --> viral and bacterial mode
if (args.mode[0] == 'bam_all'):
    print('Entering BAM  --> assembly & raw read classificaiton using Kraken and mOTUs2')
    #check 2 files per sample because they are paried end fastq files
    #if (int(args.files_per_sample[0]) != 1):
    #    raise ValueError('%s files per sample specified, this is not possible for paired end FASTQ processing' %(args.files_per_sample))
    #check name format is sample.fq.gz
    for file in inputfiles:
        if (file[-4:].lower() != '.bam'):
            raise ValueError('Files must be in sample1.bam format')
    #distribute files among working directories
    distribute_files(inputfiles, full_dir_list, 1)# args.files_per_sample[0] changed to 1 for bam files
    #navigate to each working directory and unleash ouroboros (launch seperate snakemake job launchers)
    for directory in full_dir_list:
        os.chdir(directory)
        #make scripts executable
        os.system('chmod +x all_ouroboros.py')
        #develop command to unleash ouroborus upon this world (launch fastq filtering job)
        fastq_command = ('%s "./all_ouroboros.py --outdir %s"' %(long_job, outdir))
        os.system(fastq_command)
    #revert back to original working directory
    os.chdir(cwd)




