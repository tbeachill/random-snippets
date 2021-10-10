import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import scipy

data_df = pd.read_csv('~/Physio_Agg.csv')




plot_df = data_df[['SGDName', 'AtAu', 'StSu', 'AuSu', 'AtSt']]
plot_df = plot_df.fillna(0)
mask = plot_df.isnull()
plot_df = plot_df.set_index('SGDName')

from scipy.spatial import distance
from scipy.cluster import hierarchy

correlations_array = plot_df

row_linkage = hierarchy.linkage(
    distance.pdist(correlations_array), method='average')

col_linkage = hierarchy.linkage(
    distance.pdist(correlations_array.T), method='average')

g = sns.clustermap(plot_df, row_linkage=row_linkage, col_linkage=col_linkage, method="average", figsize=(13, 13), cmap="RdBu_r", yticklabels=True)



#g = sns.clustermap(plot_df, cmap="RdBu_r", mask=mask, center=0, metric='euclidean')

a = g.dendrogram_row.reordered_ind

plt.savefig('~/test.png', figsize=(50,50))

outlist = []
for item in set(a):
    outlist.append(plot_df.index[item])
