#!/usr/bin/env python

import sys
import pandas as pd
from warnings import simplefilter # turns off an annoying warning
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
import matplotlib.pyplot as plt

# read in the output of abundanceTable.r
input_file = sys.argv[1]
df = pd.read_csv(input_file, sep="\t") 
#print(df.head(10))

# new df with first, second, and last column
df = df.drop(columns=df.columns[2:40])
#print(df.head(10))

# Pivot
df_pivot = df.pivot(index='SAMPLE_ID', columns='x1', values='RPKM')
#print(df_pivot)


# take just the top 20, and add everything else into a column called "other"
df_top20= pd.DataFrame(index=df_pivot.index)
top20=df_pivot.mean().sort_values()[-21:].index.tolist()
for gene in top20:
    df_top20[gene] = df_pivot[gene]
df_top20["Other"] = df_pivot.sum(axis=1) - df_top20.sum(axis=1)

# normalize each sample so that the counts add up to 1
df_top20_norm = df_top20.div(df_top20.sum(axis=1), axis=0) 
#print(df_top20_norm)

### visualize the data

# Transpose the DataFrame for easier plotting
#df_top20_norm_transposed = df_top20_norm.transpose()
#print(df_top20_norm_transposed)

# Plotting
plt.figure(figsize=(10, 6))
cmap = plt.cm.get_cmap('gist_rainbow', len(df_top20.columns)+1) # Specify a colormap with enough unique colors
#df_top20_norm_transposed.plot(kind='barh', stacked=True)
df_top20_norm.plot(kind='barh', stacked=True, colormap=cmap)
# plt.xlabel('Relative Abundance')
# plt.ylabel('Organisms')
plt.title('Top 20 Genes RPKM')
plt.legend(title='Sample', loc='lower center', ncols=2, bbox_to_anchor=(0.5, -0.8))

# Save as PDF
plt.savefig('genes_rpkm.pdf', bbox_inches='tight')