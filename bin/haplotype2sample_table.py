#!/usr/bin/env python3

# Created on Ap 02 2021
# @author: alineBini

import pandas as pd
import re
import csv
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input_files", required=True, help= "Input files", nargs='+')
ap.add_argument("-p", "--primer", required=True, help= "Primer name")
ap.add_argument("-c", "--min_copies", required=True, help= "default minimum copy of reads per haplotype")
args = vars(ap.parse_args())

all_files = args["input_files"]
primer = args["primer"]
min_copies = args["min_copies"]
min_copies = int(min_copies)

r = re.compile(".+_%s" % primer)
files = list(filter(r.match, all_files))

haplotype = pd.DataFrame()
i = 0

for file in files:

    sample = re.findall('(.*)_%s' % primer, file, re.IGNORECASE)[0]
    collapsed = pd.read_csv(file, sep='\n', header=None, lineterminator='>', names=['Haplotype', 'seq', 'rem'])
    
    if(collapsed.size == 0):
        continue
    else:
        collapsed.drop(columns=['rem'], inplace=True)
        new = collapsed['Haplotype'].str.split("-", expand=True)
        collapsed['Haplotype'] = new[0]
        collapsed[sample] = new[1].astype(int)

        if i == 0:
            haplotype = collapsed
        else:
            haplotype = pd.merge(haplotype, collapsed[['seq', sample]], on='seq', how = 'outer')

        haplotype.fillna(0, inplace=True)        
        i = i + 1

if(haplotype.size != 0):
    haplotype['total'] = haplotype.sum(axis=1)
    haplotype = haplotype[(haplotype['total'] > min_copies)]
    haplotype.sort_values(by=['total'], ascending=False, inplace=True)
    haplotype.iloc[:, 2:] = haplotype.iloc[:, 2:].astype(int)
    haplotype['Haplotype'] = range(1, haplotype.shape[0]+1)
    haplotype.drop(columns=['seq', 'total'], inplace=True)
    haplotype.to_csv(primer+'.haplotype2sample.txt', sep='\t', index=False, header=True, quoting=csv.QUOTE_NONE)
