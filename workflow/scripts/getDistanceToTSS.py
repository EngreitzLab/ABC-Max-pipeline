import pandas as pd
import numpy as np
import sys 

data = pd.read_csv(sys.argv[1], sep="\t")

genes = pd.read_csv("/users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.bed", sep="\t", header=None)
genes['midpoint'] = genes[1]+250

#if 'DistanceToTSS' in list(data.columns) or 'distance' in list(data.columns):
#    data[['DistanceToTSS']].to_csv(sys.argv[3], sep="\t", index=False, header=False) 
#else:
data['TargetGeneTSS'] = 0
data = data.loc[data['TargetGene'].isin(list(genes[3]))]
data = data.dropna()
data = data.sample(frac=1)
data = data.sample(n=100000)
genes = genes.loc[genes[3].isin(list(data['TargetGene']))]
for gene in genes[3].drop_duplicates():
    matches = genes.loc[genes[3]==gene].index.astype('int')
    data_matched = data.loc[data['TargetGene']==gene].index.astype('int')
    data.loc[data_matched, 'TargetGeneTSS'] = genes.loc[matches[0], 'midpoint']

data['midpoint'] = data['start']+0.5*(data['end']-data['start'])
data['midpoint'] = data['midpoint'].astype('int')
data['distanceToTSS'] = np.abs(data['midpoint'] - data['TargetGeneTSS'])
data.to_csv(sys.argv[2], sep="\t", index=False)
#    data[['distanceToTSS']].to_csv(sys.argv[3], sep="\t", index=False, header=False)
