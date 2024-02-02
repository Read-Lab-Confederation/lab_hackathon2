#!/usr/bin/env python3
import sys
import pandas as pd

# Get the file from terminal
tsvfile = sys.argv[1]

# Read the tsv file into a dataframe
data=pd.read_csv(tsvfile,sep='\t')

# Count the number of rows (lines)
num_rows = len(data)

# Count the number of columns
num_columns = len(data.columns)

# Print the results
print(f'Number of rows: {num_rows}')
print(f'Number of columns: {num_columns}')