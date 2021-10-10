import pandas as pd
import numpy as np

data = pd.read_csv('~/2017_Aggregate_Assay/2019_MBR_All_Sep/proteinGroups.csv')

to_average = []
inagg_list = []

# Take the mean of the highest 2 peptide counts and append to the aggregate list
i=0
while i < len(data['Unique peptides 04']):
    num1 = sorted([int(data['Unique peptides 05'][i]), int(data['Unique peptides 11'][i]), int(data['Unique peptides 17'][i])])[2]
    num2 = sorted([int(data['Unique peptides 05'][i]), int(data['Unique peptides 11'][i]), int(data['Unique peptides 17'][i])])[1]
    if np.mean([num1, num2]) > 2:
        inagg_list.append(data['Protein IDs'][i].split('_')[0])
    i += 1
    
# Create a list with all detected protein names    
insol = []    
for item in data['Protein IDs']:
    insol.append(item.split('_')[0])

# Remove all names that are present in the agg list from the sol list
n=0
while n < 5:
    i=0
    for item in insol:
        if item in inagg_list:
            insol.pop(insol.index(item))
            i = i - 1
        i += 1
    n += 1
    
# Now want to attach whether or not a protein is present in aggregates with
# The physiochemical properties dataset.

physio = pd.read_csv('~/All_Sep_Physio.csv')

# Create a dictionary of whether or not a protein is aggregated or not
agg_dict = {}
for item in insol:
    agg_dict[item] = 'No'
    
for item in inagg_list:
    agg_dict[item] = 'Yes'
    
physio['Aggregated'] = ''

i=0
for item in data['Protein IDs']:
    physio['Aggregated'][i] = agg_dict[item.split('_')[0]]
    i += 1
    

physio.to_csv('~/Physio_Agg.csv')