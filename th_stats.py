#!/usr/bin/env python

"""
Author: Geert van Geest
"""

from __future__ import division #if python2
import sys
# convert tophat statistics to a dic with chromosomenumbers as keys and lists of 
# reads per kilobase per sample as values
def conv_tab(statfile):
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
  return statdic
   
if __name__ == "__main__":
  #define statistics file
  wd = sys.argv[1]
  statfile = wd + '/txdout/tophat/tophat_statistics.txt'
  statdic = conv_tab(statfile)
  #write table in tophat output as conv_stat_table.txt
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

  '''
  COMMAND LINE = python /local/data/course/project/groups/mango/project_final/th_stats.py /local/data/course/project/groups/mango/project_final
  '''