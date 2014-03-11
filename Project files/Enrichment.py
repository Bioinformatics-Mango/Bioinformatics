#!/usr/bin/env python
"""
Name: Nick Alberts
Function: retrieve ensembl_transcript_id for enrichment
"""
import sys
#Retrieve ensmbl list from diff file extension
ensembl_id=[]
with open((sys.argv[1])) as infile:
  content=infile.readlines()
  infile.close
for lines in content[1:]:
  Line=lines.split()
#Selects values between 0 and 0.05. Tab 11 contains p-values;Tab 12 contains q-values  
  if 0 <= float(Line[11]) <= 0.05:
        if 0 <= float(Line[12]) <= 0.05:
          ensembl_id.append(Line[2])
  del Line

#ensembl_id to new file for further processing
with open((sys.argv[2]),'w') as outfile:
  for transcript_ids in ensembl_id:
    outfile.write(transcript_ids+'\n')
  outfile.close

