#!usr/bin/env python

'''Author: Nikos
Script for running tophat for files in a directory'''

'''Authors: Nikos, Marcel
Script that reads the filenames in a given directory
and returns a list of tuples with related couples of files
'''

#IMPORTS#
import re
import os
import sys
import subprocess

#FUNCTIONS
def make_file_pairs(path):
  '''Takes the sequence reads files in fastq format from a directory
  that contains them (format: conditionx_1.fastq,conditionx_2.fastq for each condition)
  and returns a list of tuples containing the related pairs
  [(condition1_1, condition1_1), ...,(conditionx_1,conditionx_2)]
  '''
  file_list = os.listdir(path)
  
  #Matches only files starting with paired
  regForward = re.compile('^(paired_)(.*_1\.fastq)') #Subject to change
  regReverse = re.compile('^(paired_)(.*_2\.fastq)') #Subject to change
  
  #Create two lists: 1 for forward and 1 or the reverse strand reads#
  lForward = []
  lReverse = []
  for filename in file_list:
    if (len(regForward.findall(filename))>=1):
      lForward.append(filename)
    elif (len(regReverse.findall(filename))>=1):
      lReverse.append(filename)
      
  #Sorts the lists to make 1:1 matches for the pairing    
  lForward.sort()
  lReverse.sort()
  
  #Checks the lists are of the same length and creates the pairs.
  #Else, raises an error
  pairs_list = []
  
  if len(lForward) != len(lReverse):
    raise ValueError('Please provide a list of sequences of paired-end reads')
  else:
    for i in range(len(lForward)):
      pairs_list += [(lForward[i], lReverse[i])]
  return pairs_list

def run_tophat2(reads_folder):
  genes = path + 'Saccharomyces_cerevisiae/Ensembl/EF2/Annotation/Genes/genes.gtf'
  bowtie2index = path + 'Saccharomyces_cerevisiae/Ensembl/EF2/Sequence/Bowtie2Index/genome'
  pairs = make_file_pairs(reads_folder)
  for pair in pairs:
    file1 = pair[0]
    file2 = pair[1]
    output_folder = path + 'project/tophat_testrun/'+ file1[7:-8]
    cmd = 'tophat2 -G '+ genes+' '+ '-o '+output_folder + ' ' + bowtie2index + ' '+  reads_folder+file1 + ' ' + reads_folder+file2
    print cmd+'\n'
    #subprocess.check_call(cmd, shell = True)
    
    
# reading files in directory, creating pairs while running
# list [match1,match2]
#  for each match create output dir with first part of file name
#  run tophat output is created dir
#  create 2 subfolers in output directory, 1 for pipeline 1 for manual usage
#  move files to right folders

    
###MAIN###
if __name__ == '__main__':
  path = '/local/data/course/project/groups/mango/'
  reads_folder = path + 'project/improved_quality_reads_crop15/'
  run_tophat2(reads_folder)