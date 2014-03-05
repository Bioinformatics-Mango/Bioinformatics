#!/usr/bin/env python
"""
Name: Nick Alberts
Function: Gene id for enrichment
"""
import sys
#retrieve ensmbl list from gff or gtf DONE
ensembl_id=[]
with open((sys.argv[1])) as infile:
  content=infile.readlines()
  infile.close
for lines in content[1:]:
  Line=lines.split()
  if Line[13] == 'yes':
#if (pvalue bottom) <= Line[11]  <= (pvalue upper)
    ensembl_id.append(Line[2])
  del Line
ensembl_id

#ensembl_id to new file for further processing
with open((sys.argv[2]),'w') as outfile:
  for ids in ensembl_id:
    outfile.write(ids+'\n')
  outfile.close

#input ensmbl list in david
command=['wget http://david.abcc.ncifcrf.gov/api.jsp?type=ENSEMBL_GENE_ID&ids=','&tool=geneReportFull']
data=[]

for ids in ensembl_id[:-1]:
  data.append(ids.replace(ids,ids+','))
x=''.join(data)
command.insert(1,ensembl_id[-1])
command.insert(1,x)
web_api=''.join(command)
print web_api

#subprocess.call(web_api)



#pick DAVID options
#wget http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids=2919,6347,6348,6364&tool=geneReportFull
#type ENSEMBL_GENE_ID
#annot  = a list of desired annotation  categories separated by ","
#ids  = a list of user's gene IDs separated by ","
#tool  = one of DAVID tool names  summary
