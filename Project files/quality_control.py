
#!/usr/bin/env python

"""
Author: Geert van Geest
This script takes a workign directory, a .txt file with new filenames and a path to the folder with the fastq files
It creates three new folders containing:
- fastqc output of the raw sequences (raw_fastqc)
- trimmed fastq files (improved_quality_reads)
- fastqc output of trimmed files (trimmed_fastqc)
"""

from __future__ import division #if python2
import sys
from subprocess import check_output
import os

'''Takes a directory and return as alist with all fastq files including their paths in that directory including subdirectories
'''
def list_fastq_files(directory):
  fils = {}
  for (dirpath, dirnames, filenames) in os.walk(directory):
    fils[dirpath] = filenames  
  for k, v in fils.items():
    path_files = []
    for it in v:
      path_files.append(k+'/'+it)
  fastq_files = []
  for it in path_files:
    if it.endswith('.fastq'):
      fastq_files.append(it)
  return fastq_files

  """Takes a .txt file with new filenames in alphabetical order of the old filenames and returns a list of the names
  """
def list_new_filenames(new_names_file):
  f = open(new_names_file, 'rU')
  list_new = []
  for line in f.readlines():
    list_new.append(line.strip())
  return list_new
  
'''Takes the path of the working directory and the name of a new directory made in the wd. It makes this path if it doesn't exist yet
It returns the path of the new directory
'''
def make_dir(wd, dirname):
  newdir = wd + dirname
  if not os.path.exists(newdir):
    os.makedirs(newdir)
  return wd + dirname

  '''Takes a list of paths to files and the directory for the output
  Runs fastqc on the files and puts them in the output directory
  '''
def run_fastqc(dirlist, outdir):
  fastqc_cmd = ['fastqc'] + dirlist + ['-o', outdir]
  joined_cmd = ' '.join(fastqc_cmd)
  runfastqc = check_output(joined_cmd, shell=True) 
  
  '''Takes a path to a file, output directory and output filenames
  runs fastq_quality_trimmer on these files
  '''
def run_trimmer(infile, outdir, outfile):
  outdirfile = outdir + outfile
  filter_cmd_list = ['fastq_quality_trimmer -t 30 -l 25 -i' , infile,'-o', outdirfile]
  filter_cmd = ' '.join(filter_cmd_list)
  filt = check_output(filter_cmd, shell = True)

if __name__ == "__main__":
  wd = sys.argv[1] #specify workign directory
  new_filenames = sys.argv[2] #.txt file of new filenames
  seqpath = sys.argv[3] #specify the path to the folder containing the raw sequences '/local/data/course/project/RNAseq/yeast/'
  new_name_list = list_new_filenames(new_filenames)
  
  # make a list of paths to input files
  raw_data_paths = list_fastq_files(seqpath)
  
  # run the fastqc command on all raw files in directory
  qc_ini = make_dir(wd , 'raw_fastqc')
  run_fastqc(raw_data_paths, qc_ini)
  
  # trim all raw fastq files
  trimdir = make_dir(wd, 'improved_quality_reads')
  for i in range(len(raw_data_paths)):
    run_trimmer(raw_data_paths[i], trimdir,  new_name_list[i])
  
  # run fastqc on all trimmed files
  qc_trimmed = make_dir( wd, 'trimmed_fastqc')
  trimmed_reads_paths = list_fastq_files(trimdir)
  run_fastqc(trimmed_reads_paths, qc_trimmed)
  
  '''
  COMMANDLINE: python Quality_control.py /local/data/course/project/groups/mango/project/ new_filenames.txt /local/data/course/project/RNAseq/yeast/
  '''

