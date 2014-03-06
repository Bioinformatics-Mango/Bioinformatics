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

import re
import numpy

#import urllib2
## end imports
##################################################################

## Functions
def callCuffDiff():
    ''' Calls the cuffdiff program, if not able to run the program it will give an error message
    '''
    try:    
        if(len(subprocess.check_output(["which","cuffdiff"]))>0):
            print('cuffdiff is installed')
            bPermission = checkPermissions('test.txt', ['R','W','X'])
            if(not bPermission):
                print(bPermission)
            
            bPermission = checkPermissions('test', ['R','W','X'])
            if(not bPermission):
                print(bPermission)
            bPermission = checkPermissions('test2', ['R','W','X'])
            if(not bPermission):
                print(bPermission)
            bPermission = checkPermissions('1', ['R','W','X'])
            if(not bPermission):
                print(bPermission)
            bPermission = checkPermissions('1/a.txt', ['R','W','X'])
            if(not bPermission):
                print(bPermission)
            bPermission = checkPermissions('3', ['R','W','X'])
            if(not bPermission):
                print(bPermission)
    	else:
    	  print('cuffdiff is not installed')
    except OSError as e:
        print(e)
    return ''


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
    
    callCuffDiff()
    
    
        

## End Console call    