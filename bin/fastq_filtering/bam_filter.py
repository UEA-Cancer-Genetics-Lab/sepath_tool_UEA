#!/usr/bin/env python

#this script will obtain unmapped reads from host-aligned bam files

import os
import glob
import time
import argparse
import pysam
import gzip

#build arg parser here
parser = argparse.ArgumentParser(description='Obtain unmapped reads from host-aligned bam files')
parser.add_argument('--input_bam', nargs=1, required=True, help='Enter the absolute path of a bam file to process')
args = parser.parse_args()

#bamfile set
bamfile = str(args.input_bam[0])


ID=(bamfile.split('/')[-1])[:-4] #form identification number
#develop commands to obtain unmapped files
samfile = pysam.AlignmentFile(bamfile, "rb")
#set prefix without any extension
prefix = bamfile.split('/')[-1]
prefix = prefix.split('.bam')[0]
#Set R1 and R2 filenames
R1 = ('%s_R1.fq' %prefix)
R2 = ('%s_R2.fq' %prefix)
#Open files to writeto
outfile1 = open(R1, "w")
outfile2 = open(R2, "w")
count = 0
for read in samfile.fetch(until_eof=True):
    #save read if it is unmapped or if its mate is unmapped
    if read.is_unmapped or read.mate_is_unmapped:
        #skip read if it's not the primary alignment
        if read.is_secondary:
            count = count + 1
            if count % 100000000 == 0:
                print('Processed %s million alignments' %(count/1000000))
            continue
        if read.is_read1:
            read1 = '@' + (read.query_name) + '\n' + (read.query_sequence) +  '\n' + '+' + '\n' + ((read.qual)[:]) + '\n'
            outfile1.write(str(read1))
        if read.is_read2:
            read2 = '@' + (read.query_name) + '\n' + (read.query_sequence) +  '\n' + '+' + '\n' + ((read.qual)[:]) + '\n'
            outfile2.write(str(read2))
    count = count + 1
    if count % 100000000 == 0:
        print('Processed %s million alignments' %(count/1000000))
samfile.close()
outfile1.close()
outfile2.close()
#gzip compress fastq output files 
print('Compressing...')
os.system('gzip *.fq')


