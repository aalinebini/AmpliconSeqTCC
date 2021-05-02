#!/usr/bin/env python3

# Created on Ap 01 2021
# @author: alineBini

import pandas as pd
import csv
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input_file", required=True, help= "Path of the input file")
args = vars(ap.parse_args())

input = args["input_file"]

file = pd.read_csv(input, sep='\n', header=None, lineterminator='>', names=['id', 'seq', 'rem'])
file.drop(columns=['rem'], inplace=True)
new = file['id'].str.split("-", expand=True)
file['id'] = new[0]
file['copies'] = new[1]
file = file[['id', 'copies', 'seq']]

file.to_csv(input+'.haplotypes.txt', sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)
