#!/usr/bin/env python

#this script will returnm a list of files with given extensions from a directory as input_files.txt

import glob
import os

cur_dir = (os.getcwd() + "/")
finalfile = ("%sinput_files.txt" %cur_dir)

if os.path.isfile((cur_dir) + "input_files.txt") == True:
    del_file = input("There is already an input_files.txt in the current directory, overwrite? (y/n) ")
    del_file = del_file.lower()
    if del_file == "y":
        os.system("rm -f input_files.txt")
    else:
        print("no selected, exiting script")
        raise Exception("input_files.txt already present in working directory, exiting script")    

#return absolute path of input directory
dir = os.path.abspath(input("Enter directory containing files: "))

ext = input("Enter extension of files to process (i.e. .fastq.gz): ")

searchterm = "%s/*%s" % (dir, ext)

FILES = glob.glob(searchterm)

if len(FILES) >=1:
    print("%s files found, writing to input_files.txt" %len(FILES))
elif len(FILES) == 0:
    raise Exception("No files found matching that criteria in that directory") 

#print files to a file in the current working directory
with open(finalfile, 'w') as writefile:
    for file in FILES:
        writefile.write(file + '\n')
    writefile.close


#put closely matching file names next to one another by sorting alphanumerically
os.system("sort input_files.txt > temp_file_list; mv temp_file_list input_files.txt")
