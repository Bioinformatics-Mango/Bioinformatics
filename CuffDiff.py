#!/usr/bin/env python
"""
Name: Marcel van de Streek
Script: CuffDiff call
"""

##################################################################
## Imports
from __future__ import division
from __future__ import print_function
#from Bio import SeqIO
#import Bio.ExPASy
#from __future__ import with_statement
from sys import argv
import operator
import string
import subprocess
from subprocess import Popen, PIPE
import os
import sys

import re
import numpy

#import urllib2
## end imports
##################################################################

## Functions
def callCuffDiff(sPath):
    ''' Calls the cuffdiff program, if not able to run the program it will give an error message
    '''
    try:
      if(checkRequiredFolders(sPath)): #All required files and folders are ready, creating cuffdiff commandline
            exit
            sCuffDiffCommanLine = 'cuffdiff'
            #Output directory
            sCuffDiffCommanLine += ' -o ' + sPath + 'txdout/cuffdiff'
            #Reference genome
            sCuffDiffCommanLine += ' -b' + sPath + 'ref_genome/genome.fa'
            #merged.gtf file  
            sCuffDiffCommanLine += ' -u ' + sPath + 'txdout/cuffmerge/merged_asm/merged.gtf'
            #batch bam directory
            #Create list of bam files for batch
            sCuffDiffCommanLine += ' '
            
            regTophadOut = re.compile('(.*)\d_thout')
            dictLabels = {}
            for folder in os.listdir(sPath+'txdout/tophat'):
                if(len(regTophadOut.findall(folder))>=1):
                    oLabelMatch = regTophadOut.match(folder)
                    dictLabels[oLabelMatch.group(1)] = ''

            #Labels
            sCuffDiffCommanLine += ' -L '
            for k,v in dictLabels.items(): #Happens twice
                sCuffDiffCommanLine += k+','
            #remove trailing komma - add space
            sCuffDiffCommanLine = sCuffDiffCommanLine[:-1] + ' '
                    
            for k,v in dictLabels.items(): #Happens twice
                regOutFile = re.compile('('+k+'.*)_thout')
                for folder in os.listdir(sPath+'txdout/tophat/'):

                    if(len(regOutFile.findall(folder))>=1):
                        oFileMatch = regOutFile.match(folder)             
                        oFileMatch = oFileMatch.group(1)
                        sCuffDiffCommanLine += sPath + 'txdout/tophat/' + folder + '/' + 'accepted_hits_'+oFileMatch+'.bam,'
                #remove trailing komma
                sCuffDiffCommanLine = sCuffDiffCommanLine[:-1] + ' '                          
            
            
#             regTophadOut = re.compile('(.)*_thout')
#             dictLabels = {}
#             for folder in os.listdir(sPath+'txdout/tophat'):
#                 #Start grouping them based on similar names, those will become the labels
#                 
#                 if(len(regTophadOut.findall(folder))>=1):
#                     sCuffDiffCommanLine += sPath + 'txdout/tophat/' + folder + '/' + 'accepted_hits.bam'
# 
#             #remove trailing komma
#             sCuffDiffCommanLine = sCuffDiffCommanLine[:-1]
#             
#             #Create list of bam files for chemostat
#             sCuffDiffCommanLine += ' '
#             
#             
#             regTophadOut = re.compile('.*_thout')
#             for folder in os.listdir(sPath+'txdout/chemostat'):
#                 if(len(regTophadOut.findall(folder))>=1):
#                     sCuffDiffCommanLine += sPath+'txdout/chemostat'+'/'+folder+'/accepted_hits.bam,'
#             #remove trailing komma
#             sCuffDiffCommanLine = sCuffDiffCommanLine[:-1]
            
            
            
            print(sCuffDiffCommanLine)
            sCuffDiffCommanLine = sCuffDiffCommanLine.strip()
            os.system(sCuffDiffCommanLine)
            #Create list of bam files for chemostat
            
            
    except OSError as e:
        print(e)
    return ''

def checkRequiredFolders(sPath):
    if(len(subprocess.check_output(["which","cuffdiff"]))>0):
        #Checking required imput files
        sErrorMessage = ''
        bPermission = checkPermissions(sPath + 'txdout/cuffdiff', ['R','W','X'])
        if(not bPermission):
            sErrorMessage += bPermission
        bPermission = checkPermissions(sPath + 'ref_genome/genome.fa', ['R'])
        if(not bPermission):
            sErrorMessage += bPermission
        bPermission = checkPermissions(sPath + 'txdout/cuffmerge/merged_asm/merged.gtf', ['R'])
        if(not bPermission):
            sErrorMessage += bPermission
            
        #Check if the bam folders are present
        # Loop through Tophat folder. get the _thout folders
        bPermission = checkPermissions(sPath + 'txdout/tophat', ['R'])
        if(not bPermission):
            sErrorMessage += bPermission
        else:
            regTophadOut = re.compile('(.*)_thout')
            for folder in os.listdir(sPath+'txdout/tophat'):
                if(len(regTophadOut.findall(folder))>=1):
                    oLabelMatch = regTophadOut.match(folder)
                    oLabelMatch = oLabelMatch.group(1)
                    bPermission = checkPermissions(sPath + 'txdout/tophat/' + folder + '/' + 'accepted_hits_'+oLabelMatch+'.bam', ['R'])
                    if(not bPermission):
                        sErrorMessage += bPermission
                    bPermission = checkPermissions(sPath + 'txdout/tophat/' + folder + '/' + 'accepted_hits_'+oLabelMatch+'.bam.bai', ['R'])
                    if(not bPermission):
                        sErrorMessage += bPermission                    
        if(not len(sErrorMessage)==0): 
            print(sErrorMessage)
            sys.exit()
    else:
        print('cuffdiff is not installed')   
    return True

def checkPermissions(sPath, lPermissions):
    """ checkPermissions checks if this script can find the file/folder and has the correct permissions
        When no problems occur, TRUE is returned. Otherwise a descriptive error message is returned
        sPath         -> Absolute path to file or directory
        lPermissions  -> List of permissions to check possible options (R,W,X)
                         R = Read, W = Write, X = Executable
        Return        -> True if file is accessible with the requested permissions
                         String with descriptive error message when it is not accessible with the requested permissions
    """ 
    if(os.path.exists(sPath)):
        #Check file permissions
        sErrorMessage = ''
        for sPermissiosn in lPermissions:
            if(sPermissiosn =='R' and not os.access(sPath, os.R_OK)):
                sErrorMessage = 'Could not open ' + sPath + ' for reading. \n'
            elif(not os.access(sPath, os.W_OK)):
                sErrorMessage += 'Could not open ' + sPath + ' for writing. \n'
            elif(not os.access(sPath, os.X_OK)):
                sErrorMessage += 'Could not open ' + sPath + ' for executing.'
        if(len(sErrorMessage)==0):
            
            return True
        else:
            return sErrorMessage
    else:
        print('Could not find: ' + sPath)
        print('Are you missing the absolute path?')
        return 'Could not find: ' + sPath   
    return ''


## End Functions

## Console call
if __name__ == "__main__":
    
    callCuffDiff('/local/data/course/project/groups/mango/project_final/')
    
    
        

## End Console call    
