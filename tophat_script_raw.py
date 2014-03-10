#!usr/bin/env python


'''Authors: Nikos, Marcel, Adithi
Script that reads the filenames in a given directory
and returns a list of tuples with related couples of files.
The list of read pairs is then passed into tophat for mapping.
The remaining unpaired reads from the trimming process are also mapped and
merged with the resulting paired reads.
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
  UregForward = re.compile('^(unpaired_)(.*_1\.fastq)') #Subject to change
  UregReverse = re.compile('^(unpaired_)(.*_2\.fastq)') #Subject to change
  
  #Create two lists: 1 for forward and 1 or the reverse strand reads#
  lForward = []
  lReverse = []
  ulForward = []
  ulReverse = []
  
  for filename in file_list:
    if (len(regForward.findall(filename))>=1):
      lForward.append(filename)
    elif (len(regReverse.findall(filename))>=1):
      lReverse.append(filename)
    elif (len(UregForward.findall(filename))>=1):
      ulForward.append(filename)
    elif (len(UregReverse.findall(filename))>=1):
      ulReverse.append(filename)
      
  #Sorts the lists to make 1:1 matches for the pairing    
  lForward.sort()
  lReverse.sort()
  ulForward.sort()
  ulReverse.sort()
  
  
  #Checks the lists are of the same length and creates the pairs.
  #Else, raises an error
  pairs_list = []
  
  if len(lForward) != len(lReverse):
    raise ValueError('Please provide a list of sequences of paired-end reads')
  else:
    for i in range(len(lForward)):
      pairs_list += [(lForward[i], lReverse[i])]
  return pairs_list,ulForward,ulReverse

def run_tophat2(reads_folder):
  genes = path + 'Saccharomyces_cerevisiae/Ensembl/EF2/Annotation/Genes/genes.gtf'
  bowtie2index = path + 'Saccharomyces_cerevisiae/Ensembl/EF2/Sequence/Bowtie2Index/genome'
  pairs,ulForward,ulReverse = make_file_pairs(reads_folder)
  #print ulForward,ulReverse 
  i = 0
  
  for pair in pairs:
    file1 = pair[0]
    file2 = pair[1]
    
    output_folder = path + 'project_final/txdout/tophat/'+ file1[7:-8]+'_thout'
    cmd = 'tophat2 -p 8 -G '+ genes+' '+ '-o '+output_folder + ' ' + bowtie2index + ' '+  reads_folder+file1 + ' ' + reads_folder+file2
    
    #print cmd+'\n'
    
    subprocess.check_call(cmd, shell = True)
    
    unpaired1 = str(ulForward[i])
    unpaired2 = str(ulReverse[i])
  
    #print unpaired1, unpaired2
  
    output_folder_unpaired1 = output_folder + '/'+unpaired1[:-6]
    cmd_unpaired1 = 'tophat2 -p 8 -G '+ genes+' '+ '-o '+output_folder_unpaired1 + ' ' + bowtie2index + ' '+  reads_folder+unpaired1
    
    output_folder_unpaired2 = output_folder +'/'+unpaired2[:-6]
    cmd_unpaired2 = 'tophat2 -p 8 -G '+ genes+' '+ '-o '+output_folder_unpaired2 + ' ' + bowtie2index + ' '+  reads_folder+unpaired2
    
    cmd_merge = 'samtools merge '+ output_folder + '/accepted_hits_' + file1[7:-8] + '.bam ' + output_folder + '/accepted_hits.bam ' + output_folder_unpaired1 + '/accepted_hits.bam ' + output_folder_unpaired2 + '/accepted_hits.bam'
    
    i += 1   
    
    #print cmd_unpaired1+'\n'
    #print cmd_unpaired2+'\n'
    #print cmd_merge+'\n'
    
    
    subprocess.check_call(cmd_unpaired1, shell = True)
    subprocess.check_call(cmd_unpaired2, shell = True)
    subprocess.check_call(cmd_merge, shell = True)
    
    
###MAIN###
if __name__ == '__main__':
  path = '/local/data/course/project/groups/mango/'
  reads_folder = path + 'project_final/trimmed_fastq/'
  run_tophat2(reads_folder)
