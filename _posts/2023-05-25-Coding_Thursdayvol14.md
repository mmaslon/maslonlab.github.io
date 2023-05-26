---
title: Processing rMATS output 
date: 2023-05-18
categories: ["Coding Thursday"]
background: /assets/theme/images//as.png
---

rMATS outputs several files, including txt files for all considered splicing events (i.e. SE (skipped exon), MXE (mutually exclusive exons), A3SS (alternative 3' splice site), A5SS (alternative 5' splice site), RI (retained intron)). There are two versions of these files, with exentsions JC or JCEC, which include only reads that span junctions or both reads that span junctions and reads that do not cross an exon boundary; respectively.

For downstream analysis, I would like to combine these files into one file. I will create two files: one filteres on FDR, one not. I will also select columns that I find useful. I will use python script to do so:

```python
#!/bin/bash/usr
import glob
import pandas as pd
import numpy as np
import os
from pathlib import Path

AS_events = 'path/to/txtfiles/with/JCevents'

#create lists for all events of FDR filtered events
events_list_all = []
events_list_sig = []

all_AS_events=sorted(glob.glob(AS_events + '/*.txt'))

for a in all_AS_events:
    df = pd.read_table(a, sep='\t', engine='python')
    #get event name from the filename, by extracting substring before last period (.txt) and splitting it further, at it as a new column to a dataframe  
    samplename_working=Path(a).stem
    samplename_working2=os.path.splitext(samplename_working)[0]
    samplename_final=os.path.splitext(samplename_working2)[0]
    df['Group']=samplename_final
    
    #subset only required columns
    df=df[['Group','GeneID','chr','strand','FDR','IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2','IncLevel1','IncLevel2','IncLevelDifference']]
    #subset only significant events
    df_sigw=df[df['FDR']<= 0.05]
    #append to appropriate lists
    events_list_all.append(df)
    events_list_sig.append(df_sigw)

#combine the lists
WT_all = pd.concat(events_list_all)
WT_sig = pd.concat(events_list_sig)

#replace NA with 0
WT_all = WT_all.fillna(0)
WT_sig = WT_sig.fillna(0)

#specify file names for these new files
sig = 'events_sig.csv'
all = 'events_all.csv'
output_sig = os.path.join(AS_events, sig)
output_all = os.path.join(AS_events, all)

#save
WT_sig.to_csv(output_sig, sep=',', header=True)
WT_all.to_csv(output_all, sep=',', header=True)

```
