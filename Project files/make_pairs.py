#!/usr/bin/env python

'''Author: Nikos, Marcel
Script that reads the filenames in a given directory
and returns a list of tuples with related couples of files
'''

#IMPORTS#
import re
import os
import sys

#FUNCTIONS
def make_file_pairs(path):
  '''Takes the sequence reads files in fastq format from a directory
  that contains them (format: conditionx_1.fastq,conditionx_2.fastq for each condition)
  and returns a list of tuples containing the related pairs
  [(condition1_1, condition1_1), ...,(conditionx_1,conditionx_2)]
  '''
  file_list = os.listdir(path)

  regForward = re.compile('(.*_1\.fastq)') #Subject to change
  regReverse = re.compile('(.*_2\.fastq)') #Subject to change
  
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
  
#MAIN#   
if __name__ == '__main__':
  #Just an example#
  path = raw_input('Please provide a directory with paired-end sequence reads files: ')
  print make_file_pairs(path)