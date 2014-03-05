
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
  path_files = []
  for k, v in fils.items():  
    for it in v:
      path_files.append(k+'/'+it)
  fastq_files = []
  for it in path_files:
    if it.endswith('.fastq'):
      fastq_files.append(it)
  return sorted(fastq_files)

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
def run_fastqc(fastq_file, outdir):
  fastqc_cmd = ['fastqc', fastq_file,'-o', outdir]
  joined_cmd = ' '.join(fastqc_cmd)
  runfastqc = check_output(joined_cmd, shell=True) 

  '''Takes the working directory and downloads and unzips Trimmomatic-0.32
  '''
def install_Trimmomatic(wd):
  install_cmd_list = ['wget www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip -O',wd+'Trimmomatic-0.32.zip']
  install_cmd = ' '.join(install_cmd_list)
  run_install = check_output(install_cmd, shell=True)
  unzip_cmd = ' '.join(['unzip', wd + 'Trimmomatic-0.32.zip'])
  run_unzip = check_output(unzip_cmd, shell=True)
  
  '''Takes paths to forward and reverse file, trimmomatic app, fasta file with adapters, output directory output filenames
  runs Trimmomatic
  '''
def run_Trimmomatic(infile_1, infile_2, path_to_trimmomatic, path_to_adapters, outdir, outfile_1, outfile_2):

  Tri_start = ['java -jar', path_to_trimmomatic, 'PE']
  Tri_files = [infile_1, infile_2, outdir + 'paired_' + outfile_1, outdir + 'unpaired_'+ outfile_1, outdir + 'paired_' + outfile_2, outdir + 'unpaired_'+outfile_2]
  Tri_options = ['ILLUMINACLIP:'+path_to_adapters+':2:30:10','HEADCROP:15', 'LEADING:30', 'TRAILING:30', 'MINLEN:25'] 
  Trimmomatic_cmd_list = Tri_start + Tri_files + Tri_options  
  Trimmomatic_cmd = ' '.join(Trimmomatic_cmd_list)
  Trim = check_output(Trimmomatic_cmd, shell=True)

if __name__ == "__main__":
  wd = sys.argv[1] #specify workign directory
  new_filenames = sys.argv[2] #.txt file of new filenames
  seqpath = sys.argv[3] #specify the path to the folder containing the raw sequences 
  new_name_list = list_new_filenames(new_filenames)
  # make a list of paths to input files
  raw_data_paths = list_fastq_files(seqpath)
  
  #install Trimmomatic if not installed yet
  Trimm_dir = wd + 'Trimmomatic-0.32'
  if not os.path.exists(Trimm_dir):
    install_Trimmomatic(wd)
        
  # run the fastqc command on all raw files in directory
  qc_ini = make_dir(wd , 'raw_fastqc')
  #for i in range(len(raw_data_paths)):
    #run_fastqc(raw_data_paths[i], qc_ini+'/'+new_name_list[i])

  # trim all raw fastq files
  trimdir = make_dir(wd, 'improved_quality_reads')
  #for i in range(0, len(raw_data_paths), 2):
   # run_Trimmomatic(raw_data_paths[i], raw_data_paths[i+1], wd+'/Trimmomatic-0.32/trimmomatic-0.32.jar',wd+'/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa', trimdir, new_name_list[i], new_name_list[i+1])
   
  # run fastqc on all trimmed files
  qc_trimmed = make_dir( wd, 'trimmed_fastqc/')
  trimmed_reads_paths = list_fastq_files(trimdir)
  #for i in range(len(trimmed_reads_paths)):
   # run_fastqc(trimmed_reads_paths[i], qc_trimmed+'/'+new_name_list)
  
  '''
  COMMANDLINE: python Quality_control.py /local/data/course/project/groups/mango/project/ new_filenames.txt /local/data/course/project/RNAseq/yeast/
  '''
