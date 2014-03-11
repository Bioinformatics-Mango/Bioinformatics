#!/usr/bin/env python
"""
Name: Mango group
Authors:    Nick Alberts, Geert van Geest, Nikolaos Pappas, Marcel van de Streek, Adithi Varadarajan
Script:     Automated pipeline to analyze differentially expressed genes under two conditions.
            Steps in this pipeline include:
            1. Take in FASTQ sequences
            2. Quality control of the imput such as:
             - Trimming low quality sequences
             - Removing primer sequences
             - 
            3. CuffLink     - Mapping of expression data against a reference genome
            4. CuffCompare  - 
            5. CuffMerge    -
            5. CuffDiff     - Analyze for differentialy expressed genes
            6. Output a list of gene's with their corresponding data such as:
             - p-values
             - q-values (adjusted p-values according to bonferroni-hogberg)
             - log2Foldchanges
             
"""

#Imports
from __future__ import division         #if python2
from __future__ import print_function   #if python2
from subprocess import check_output
from subprocess import Popen, PIPE
import operator
import os
import re
import string
import subprocess
import sys
#End Imports

### General functions
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
        #bPermission = checkPermissions(sPath + 'txdout/cuffmerge/merged_asm/merged.gtf', ['R'])
        bPermission = checkPermissions(sPath + 'txdout/cuffcompare/cuffcmp.combined.gtf', ['R'])
        if(not bPermission):
            sErrorMessage += bPermission
            
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
### End General functions


### Quality control functions

"""
Takes:
- a .txt file with new filenames in alphabetical order of the old filenames (new_names_file)
returns a list of the names
"""
def list_new_filenames(new_names_file):
    f = open(new_names_file, 'rU')
    list_new = []
    for line in f.readlines():
        list_new.append(line.strip())
    return list_new

'''
Takes:
- a directory with fastq files (directory)
returns a list with all fastq files including their paths in that directory
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

'''
Takes: the working directory (wd)
downloads and unzips Trimmomatic-0.32
'''
def install_Trimmomatic(wd):
    install_cmd_list = ['wget www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip -O',wd+'Trimmomatic-0.32.zip']
    install_cmd = ' '.join(install_cmd_list)
    run_install = check_output(install_cmd, shell=True)
    unzip_cmd = ' '.join(['unzip', wd + 'Trimmomatic-0.32.zip'])
    run_unzip = check_output(unzip_cmd, shell=True)
 
'''
Takes:
- the working directory (wd)
- the name of a new directory (dirname).
It makes a path to the new directory if it doesn't exist yet, and returns the path of the new directory
'''    
def make_dir(wd, dirname):
    newdir = wd + dirname
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    return wd + dirname

'''
Takes:
- a list of paths to fastq files (fastq_file)
- the directory for the output (outdir)
Runs fastqc on the files and puts them in the output directory
'''
def run_fastqc(fastq_file, outdir):
    fastqc_cmd = ['fastqc', fastq_file,'-o', outdir]
    joined_cmd = ' '.join(fastqc_cmd)
    runfastqc = check_output(joined_cmd, shell=True) 

'''
Takes:
- paths to forward and reverse file (infile_1; infile_2), 
- location of trimmomatic app (path_to_trimmomatic)
- fasta file with adapters (path_to_adapters)
- output directory (outdir) and output filenames (outfile_1; outfile_2)
runs Trimmomatic on infile_1 and infile_2
'''
def run_Trimmomatic(infile_1, infile_2, path_to_trimmomatic, path_to_adapters, outdir, outfile_1, outfile_2):
    Tri_start = ['java -jar', path_to_trimmomatic, 'PE']
    Tri_files = [infile_1, infile_2, outdir + 'paired_' + outfile_1, outdir + 'unpaired_'+ outfile_1, outdir + 'paired_' + outfile_2, outdir + 'unpaired_'+outfile_2]
    Tri_options = ['ILLUMINACLIP:'+path_to_adapters+':2:30:10','HEADCROP:15', 'LEADING:30', 'TRAILING:30', 'MINLEN:25'] 
    Trimmomatic_cmd_list = Tri_start + Tri_files + Tri_options  
    Trimmomatic_cmd = ' '.join(Trimmomatic_cmd_list)
    Trim = check_output(Trimmomatic_cmd, shell=True)

def qualityControlInput(workingDirectory, new_filenames,seqpath):
    try:
        new_name_list = list_new_filenames(new_filenames)
        raw_data_paths = list_fastq_files(seqpath)                              # make a list of paths to input files
        Trimm_dir = workingDirectory + 'Trimmomatic-0.32/'                      #install Trimmomatic if not installed yet
        if not os.path.exists(Trimm_dir):
            install_Trimmomatic(workingDirectory)
              
        # run the fastqc command on all raw files in directory
        qc_ini = make_dir(workingDirectory , 'fastqc_out/before_trimming/')
        for i in range(len(raw_data_paths)):
            qc_dir = make_dir(qc_ini, new_name_list[i])
            run_fastqc(raw_data_paths[i], qc_dir)
        
        # trim all raw fastq files
        trimdir = make_dir(workingDirectory, 'trimmed_fastq/')
        for i in range(0, len(raw_data_paths), 2):
            run_Trimmomatic(raw_data_paths[i], raw_data_paths[i+1], workingDirectory+'/Trimmomatic-0.32/trimmomatic-0.32.jar',workingDirectory+'/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa', trimdir, new_name_list[i]+'.fastq', new_name_list[i+1]+'.fastq')
         
        # run fastqc on all trimmed files
        qc_trimmed = make_dir( workingDirectory, 'fastqc_out/after_trimming/')
        trimmed_reads_paths = list_fastq_files(trimdir)
        trimmed_reads_filenames = []
        for it in trimmed_reads_paths:
            fil = it.split('/')[-1][:-6]
            trimmed_reads_filenames.append(fil)
          
        for i in range(len(trimmed_reads_paths)):
            qc_dir = make_dir(qc_trimmed, trimmed_reads_filenames[i])
            run_fastqc(trimmed_reads_paths[i], qc_dir)
        return True
    except:
        return sys.exc_info()[0] 
### End Quality control functions

### Tophat functions
# convert tophat statistics to a dic with chromosomenumbers as keys and lists of 
# reads per kilobase per sample as values
def conv_tab(wd):
    statfile = wd+'/txdout/tophat/tophat_statistics.txt'
    f = open(statfile)
    statdic = {}
    statdic['samplenames'] = []
    for line in f:
        row = line.split('\t')
        if row[0][0] == '@':
            statdic['samplenames'].append(line[1:].strip())
        elif row[0] != '*':
            if row[0] not in statdic.keys():
                statdic[row[0]] = [int(row[2])/int(row[1])*1000]
            else: 
                statdic[row[0]] += [int(row[2])/int(row[1])*1000]
                
    outf = open(wd + '/txdout/tophat/mapped_reads_per_1000bp.txt', 'w')
    outf.write('chr\t')
    for item in statdic['samplenames']:
        outf.write(item + '\t')
    outf.write('\n')
    for k,v in sorted(statdic.items()):
        if k != 'samplenames':
            outf.write(k +'\t')
            for item in v:
                outf.write('%s\t' % item)
        outf.write('\n')

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

def run_tophat2(path):
    #try:
        
        reads_folder = path+'trimmed_fastq/'
        genes = path + 'ref_genome/genes.gtf'
        bowtie2index = path + 'ref_genome/genome'
      #  print('start')
        pairs,ulForward,ulReverse = make_file_pairs(reads_folder)
        #print ulForward,ulReverse 
        i = 0
        
        for pair in pairs:
            file1 = pair[0]
            file2 = pair[1]
            
            output_folder = path + 'txdout/tophat/'+ file1[7:-8]+'_thout'
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
            return True                                                                         #Return true if the command is executed successfully
#     except:
#         return sys.exc_info()[0]                                                                #Return the error to the function caller if something whent wrong
#     return ''            
### End Tophat functions

### CuffDiff functions
def callCuffDiff(sPath):
    ''' Calls the cuffdiff program, if not able to run the program it will give return error message
        Return ->    True if cuffdiff was run successfully or error message when something has gone wrong
    ''' 
    try:
        if(checkRequiredFolders(sPath)): #All required files and folders are ready, creating cuffdiff commandline
            sCuffDiffCommanLine = 'cuffdiff'                                                    #Start of the cuffdiff commandline
            sCuffDiffCommanLine += ' -o ' + sPath + 'txdout/cuffdiff'                           #Absolute path to the cuffdiff output directory
            sCuffDiffCommanLine += ' -b' + sPath + 'ref_genome/genome.fa'                       #Absolute path to reference genome file
            #sCuffDiffCommanLine += ' -u ' + sPath + 'txdout/cuffmerge/merged_asm/merged.gtf'
            sCuffDiffCommanLine += ' -u ' + sPath + 'txdout/cuffcompare/cuffcmp.combined.gtf'   #Absolute path to annotation file
            sCuffDiffCommanLine += ' -p 8'                                                      #Number of CPU threads
            
            sCuffDiffCommanLine += ' --total-hits-norm'                                         #Tell Cuffdiff to count the hits per fragment, this is disabled by default
            
            #Create list of bam files and generate treatment labels based on folder names
            sCuffDiffCommanLine += ' '
            regTophadOut = re.compile('(.*)\d_thout')                                           #Regular expression that matches tophat output data per sample
            dictLabels = {}
            for folder in os.listdir(sPath+'txdout/tophat'):                                    #Loop trough list of tophat output folders
                if(len(regTophadOut.findall(folder))>=1):                                       #Only select the folders that match tophat output data for samples, others files are not relevant
                    oLabelMatch = regTophadOut.match(folder)                                    #Determine label via regular expression
                    dictLabels[oLabelMatch.group(1)] = ''                                       #Save the unique labels of the treatments in a dictionary

            sCuffDiffCommanLine += ' -L '                                                       #Adding treatment labels to the commandline
            for k,v in dictLabels.items(): #Loops twice                            
                sCuffDiffCommanLine += k+','
            sCuffDiffCommanLine = sCuffDiffCommanLine[:-1] + ' '                                #remove trailing komma - add space
                    
            for k,v in dictLabels.items(): #Loops twice                                         #Create the sample treatment groups based on the previously generated label names
                regOutFile = re.compile('('+k+'.*)_thout')                                      #Regular expression that matches tophat output data per label
                for folder in os.listdir(sPath+'txdout/tophat/'):

                    if(len(regOutFile.findall(folder))>=1):                                     #Adding each matching folder to the commandline, only 2 lists are accepted by cuffdiff
                        oFileMatch = regOutFile.match(folder)                                   #Extract foldername (1/2)
                        oFileMatch = oFileMatch.group(1)                                        #Extract foldername (2/2)
                        sCuffDiffCommanLine += sPath + 'txdout/tophat/' + folder + '/' + 'accepted_hits_'+oFileMatch+'.bam,' #Add accepted_hits bam file to the commandline
                #remove trailing komma
                sCuffDiffCommanLine = sCuffDiffCommanLine[:-1] + ' '                            #remove trailing komma - add space
            sCuffDiffCommanLine = sCuffDiffCommanLine.strip()                                   #remove trailing space
            #print(sCuffDiffCommanLine)                                                          #Showing cuffdiff commandline for debugging
            os.system(sCuffDiffCommanLine)                                                      #Execute cuffdiff commandline
            return True                                                                         #Return true if the command is executed successfully
    except:
        return sys.exc_info()[0]                                                                #Return the error to the function caller if something whent wrong
    return ''
### End CuffDiff functions

### Create Volcano plot
def createVolcanoplot(workingDirectory):
    cuffdiff_out = workingDirectory + '/txdout/cuffdiff/gene_exp.diff'
    
    
    #create directory for plots
    plotdir = make_dir(workingDirectory, '/plots')
    #execute R script 
    source_function = ["\"source('"+workingDirectory+"volcano.R')"]
    execute_function = [";volcano('" + cuffdiff_out +"','"+ plotdir +"/volcano_plot.png')\""]
    R_cmd_list = ["R -e "]+source_function+execute_function
    R_cmd = ''.join(R_cmd_list)
    check_output(R_cmd, shell=True)
    '''
    COMMAND LINE = python Volcano_plot.py '/local/data/course/project/groups/mango/project_final'
    '''
### End Volcano plot


if __name__ == "__main__":
    workingDirectory = sys.argv[1]                      #Absolutpath to the working directory for the pipeline.
    new_filenames = workingDirectory + sys.argv[2]      #'.txt file' containing a list with filenames used to rename the input sequence files. Each file is entered at a new line
    if(len(sys.argv)==4):#Optional argument
        seqpath = sys.argv[3]                           #Absolute path to the folder containing the raw sequences, default: raw_fastq
    else:
        seqpath = workingDirectory + '/raw_fastq'
    
    #bSucceeded =  qualityControlInput(workingDirectory,new_filenames,seqpath)
    #if(not bSucceeded):
    #    print(bSucceeded)
    #    exit()
  
    bSucceeded = run_tophat2(workingDirectory)
    if(not bSucceeded):
        print(bSucceeded)
        exit()   
	
    subprocess.check_call(workingDirectory+'Cuff.sh ' + workingDirectory, shell = True)    
   
    bSucceeded = callCuffDiff(workingDirectory)
    if(not bSucceeded):
        print(bSucceeded)
        exit()        
    
    createVolcanoplot(workingDirectory)
    conv_tab(workingDirectory)
    print('Congratulations! You have reached the end of the DE-pipeline. Hoepfully you have something more usefull now!')
        
                            
    
