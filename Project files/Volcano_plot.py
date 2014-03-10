#!/usr/bin/env python

"""
Author: Geert van Geest
"""

from __future__ import division #if python2
import sys
from subprocess import check_output
import os

def make_dir(wd, dirname):
  newdir = wd + dirname
  if not os.path.exists(newdir):
    os.makedirs(newdir)
  return wd + dirname

if __name__ == "__main__":
  
  #define gene_exp.diff
  wd = sys.argv[1]
  cuffdiff_out = wd + '/txdout/cuffdiff/gene_exp.diff'
  #create directory for plots
  plotdir = make_dir(wd, '/plots')
  #execute R script 
  source_function = ["\"source('/local/data/course/project/groups/mango/project_final/volcano.R')"]
  execute_function = [";volcano('" + cuffdiff_out +"','"+ plotdir +"/volcano_plot.png')\""]
  R_cmd_list = ["R -e "]+source_function+execute_function
  R_cmd = ''.join(R_cmd_list)
  check_output(R_cmd, shell=True)
  '''
  COMMAND LINE = python Volcano_plot.py '/local/data/course/project/groups/mango/project_final'
  '''